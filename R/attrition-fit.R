# Nested retention fits use full rosters; outcome learners see observed rows only.
.attrition_base <- function(d, outcome, predictors, learner, recipe_args, family) {
  if (family == "binomial") .base_fit(d, outcome, predictors, learner, recipe_args) else
    .response_base(d, outcome, predictors, learner, recipe_args, family)
}
.attrition_base_predict <- function(f, data, family) {
  if (family == "binomial") .base_predict(f, data) else .response_base_predict(f, data)
}
.attrition_prediction <- function(f, data, family, calibrated = TRUE) {
  if (family == "binomial") .pipeline_predict(f, data, calibrated) else
    .response_pipeline_predict(f, data, calibrated)
}
.attrition_calibrate <- function(y, p, w, family, method) {
  if (family == "binomial") svy_calibrate(y, p, w, method = method) else
    .response_calibrator(y, p, w, family, method)
}
.attrition_validate_fit <- function(design, outcome, predictors, learners, inner_v, family, spec) {
  z <- .attrition_roster(design, spec)
  if (!length(predictors) || !all(c(outcome, predictors) %in% names(design$data)) ||
      any(predictors %in% c(outcome, design$weights, design$ids, design$strata, design$fpc,
        design$subject, design$household, design$fold_groups, spec$response, spec$analysis)))
    stop("Predictors must exist and exclude outcome, response indicator and design identifiers.")
  if (!length(learners) || !all(vapply(learners, inherits, logical(1), "svy_learner")))
    stop("learners must be a nonempty list of svy_learner objects.")
  if (length(inner_v) != 1L || !is.finite(inner_v) || inner_v < 2L || inner_v != as.integer(inner_v))
    stop("inner_v must be an integer >= 2.")
  y <- design$data[[outcome]][z$eligible]
  if (family == "binomial") .binary(y) else .response_y(y, family)
  z
}
.attrition_pipeline <- function(d, outcome, predictors, learners, inner_v, calibration,
                                stack, recipe_args, seed, family, spec) .seed(seed, {
  z <- .attrition_roster(d, spec); rows <- which(z$eligible); n <- length(rows)
  nf <- min(inner_v, length(unique(.fold_cluster(d))))
  if (nf < 2L) stop("Too few training clusters for inner cross-fitting.")
  folds <- svy_split(d, nf, seed = seed)[, 1L]
  y <- d$data[[outcome]][rows]; lev <- levels(y)
  blank <- if (family == "multinomial") matrix(NA_real_, n, length(lev), dimnames = list(NULL, lev)) else rep(NA_real_, n)
  predictions <- rep(list(blank), length(learners)); w <- rep(NA_real_, n)
  status <- rep("ok", length(learners)); failed <- rep(FALSE, length(learners)); retention_audit <- list()
  for (k in seq_len(nf)) {
    a <- svy_attrition_weights(d, spec, training = folds != k, seed = seed + k)
    tr <- .attrition_observed(d, a$weights, z$eligible & folds != k)
    te <- which(folds[rows] == k)
    if (!length(te)) stop("An inner validation partition has no observed eligible outcomes.")
    w[te] <- a$weights[rows[te]]
    retention_audit[[k]] <- list(fold = k, audit = a$audit, training_clusters = a$model$training_clusters)
    for (j in seq_along(learners)) tryCatch({
      f <- .attrition_base(tr, outcome, predictors, learners[[j]], recipe_args, family)
      p <- .attrition_base_predict(f, d$data[rows[te], , drop = FALSE], family)
      if (family == "multinomial") predictions[[j]][te, ] <- p else predictions[[j]][te] <- p
      if (f$fit$status != "ok") status[j] <- paste(status[j], f$fit$status, sep = "; ")
    }, error = function(e) { failed[j] <<- TRUE; status[j] <<- conditionMessage(e) })
  }
  loss <- vapply(seq_along(learners), function(j) if (failed[j]) Inf else
    svy_metrics(y, predictions[[j]], w, family)[[1L]], numeric(1))
  good <- which(is.finite(loss))
  if (!length(good)) stop("All learner candidates failed inner validation: ", paste(unique(status), collapse = "; "))
  mix <- numeric(length(learners))
  if (stack && length(good) > 1L) {
    soft <- function(a) { z <- exp(c(a, 0) - max(c(a, 0))); z / sum(z) }
    opt <- stats::optim(rep(0, length(good) - 1L), function(a)
      svy_metrics(y, .response_mix(predictions[good], soft(a)), w, family)[[1L]], method = "BFGS")
    if (opt$convergence != 0) stop("Stacking optimizer did not converge.")
    mix[good] <- soft(opt$par)
  } else mix[which.min(loss)] <- 1
  cal <- .attrition_calibrate(y, .response_mix(predictions, mix), w, family, calibration)
  a <- svy_attrition_weights(d, spec, seed = seed + 10000L)
  observed <- .attrition_observed(d, a$weights, z$eligible)
  fits <- lapply(which(mix > 0), function(j) .attrition_base(observed, outcome, predictors,
    learners[[j]], recipe_args, family))
  list(fits = fits, mixture = mix[mix > 0], candidate_weights = mix, calibrator = cal,
    inner_loss = loss, inner_status = status, inner_folds = folds[rows],
    attrition_model = a$model, attrition_audit = a$audit, inner_attrition_audit = retention_audit,
    inner_attrition_weights = w, roster_inner_folds = folds)
})

.attrition_fit <- function(design, outcome, predictors, learners, v, repeats, inner_v,
                           calibration, stack, recipe_args, seed, family, spec, call) {
  z <- .attrition_validate_fit(design, outcome, predictors, learners, inner_v, family, spec)
  if (family != "binomial") calibration <- match.arg(calibration, c("auto", "none"))
  rows <- which(z$eligible); n <- length(rows); y <- design$data[[outcome]][rows]; lev <- levels(y)
  folds <- svy_split(design, v, repeats, seed); categorical <- family == "multinomial"
  raw <- pred <- if (categorical) array(NA_real_, c(n, length(lev), repeats),
    dimnames = list(NULL, lev, NULL)) else matrix(NA_real_, n, repeats)
  weights <- matrix(NA_real_, n, repeats); audit <- outer <- retention <- list()
  for (r in seq_len(repeats)) for (k in seq_len(v)) {
    tr <- which(folds[, r] != k); test_roster <- which(folds[, r] == k)
    te <- which(folds[rows, r] == k)
    if (!length(te)) stop("An outer validation partition has no observed eligible outcomes.")
    f <- .attrition_pipeline(.subset_design(design, tr), outcome, predictors, learners,
      inner_v, calibration, stack, recipe_args, seed + r * 1000L + k, family, spec)
    a <- .attrition_apply(f$attrition_model, .subset_design(design, test_roster))
    weights[te, r] <- a$weights[match(rows[te], test_roster)]
    p <- .attrition_prediction(f, design$data[rows[te], , drop = FALSE], family, FALSE)
    pp <- if (family == "binomial") stats::predict(f$calibrator, p) else .response_calibrate_predict(f$calibrator, p)
    if (categorical) { raw[te, , r] <- p; pred[te, , r] <- pp } else { raw[te, r] <- p; pred[te, r] <- pp }
    audit[[length(audit) + 1L]] <- data.frame(repeat_id = r, fold = k, train_n = sum(z$eligible[tr]),
      test_n = length(te), candidate = seq_along(learners), inner_loss = f$inner_loss,
      mixture_weight = f$candidate_weights, inner_status = f$inner_status,
      calibration_status = f$calibrator$status,
      learner_status = paste(vapply(f$fits, function(a) a$fit$status, character(1)), collapse = "; "))
    retention[[length(retention) + 1L]] <- cbind(repeat_id = r, fold = k, a$audit)
    outer[[length(outer) + 1L]] <- f
  }
  metrics <- do.call(rbind, lapply(seq_len(repeats), function(r) svy_metrics(y,
    if (categorical) matrix(pred[, , r], n, length(lev), dimnames = list(NULL, lev)) else pred[, r], weights[, r], family)))
  consensus <- if (categorical) apply(pred, c(1L, 2L), mean) else rowMeans(pred)
  final <- .attrition_pipeline(design, outcome, predictors, learners, inner_v,
    calibration, stack, recipe_args, seed + 900000L, family, spec)
  observed <- .subset_design(design, rows); observed$w <- rowMeans(weights)
  observed$data[[observed$weights]] <- observed$w; observed$repweights <- NULL; observed$replicate_df <- NULL
  structure(list(final = final, outer_fits = outer, folds = folds[rows, , drop = FALSE],
    predictions_raw = raw, predictions = pred, consensus = consensus, metrics = metrics,
    audit = do.call(rbind, audit), outcome = outcome, predictors = predictors, family = family,
    levels = lev, design = observed, diagnostics = svy_diagnostics(observed),
    attrition = list(roster = design, spec = spec, rows = rows, roster_folds = folds,
      weights = weights, audit = do.call(rbind, retention),
      inference = "Fixed-score inference also conditions on estimated retention weights; use baseline-roster replicates to refit all nuisance models",
      consensus = "Consensus predictions and mean cross-fitted weights; distinct from mean repeat metrics"),
    configuration = list(outcome = outcome, predictors = predictors, learners = learners,
      v = v, repeats = repeats, inner_v = inner_v, calibration = calibration,
      stack = stack, recipe_args = recipe_args, seed = seed, family = family, attrition = spec),
    call = call, seed = seed),
    class = c("svy_attrition_fit", if (family != "binomial") "svy_response_fit", "svy_learn_fit"))
}

#' @export
predict.svy_attrition_fit <- function(object, newdata, calibrated = TRUE, ...) {
  .attrition_prediction(object$final, newdata, object$family, calibrated)
}
#' @export
summary.svy_attrition_fit <- function(object, ...) {
  list(family = object$family, design = object$diagnostics,
    mean_repeat_metrics = colMeans(object$metrics), repeat_sd = apply(object$metrics, 2L, stats::sd),
    consensus_metrics = svy_metrics(object$design$data[[object$outcome]], object$consensus,
      object$design$w, object$family), attrition = object$attrition$audit,
    inference = object$attrition$inference, calibration_status = unique(object$audit$calibration_status))
}
