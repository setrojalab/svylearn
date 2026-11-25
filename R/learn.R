#' Specify a weighted response learner
#' @param engine glm, glmnet, ranger, xgboost, nnet or custom.
#' @param parameters Named engine arguments. glmnet defaults to alpha = 0.5 and lambda = 0.01.
#' @param fit,predict Custom functions with signatures fit(x,y,w) and predict(model,x).
#'   For multinomial fitting y is a factor with declared levels and predict must
#'   return a row-normalized probability matrix whose columns name those levels.
#'   Other responses and predictions are numeric vectors. The nnet adapter
#'   supports binomial, gaussian and multinomial; the other built-in adapters
#'   also support poisson. Multinomial glm uses baseline-category softmax with
#'   a slope ridge penalty lambda (default 1e-6), and maxit (default 5000).
#' @return A svy_learner specification. No engine is loaded at construction.
#' @export
svy_learner <- function(engine = "glm", parameters = list(), fit = NULL, predict = NULL) {
  engine <- match.arg(engine, c("glm", "glmnet", "ranger", "xgboost", "nnet", "custom"))
  if (engine == "custom" && (!is.function(fit) || !is.function(predict))) stop("Provide custom fit and predict functions.")
  structure(list(engine = engine, parameters = parameters, fit = fit, predict = predict), class = "svy_learner")
}
.fit_engine <- function(x, y, w, spec) {
  a <- spec$parameters; w <- w / mean(w)
  if (min(sum(y == 1), sum(y == 0)) < 2L) {
    neff <- sum(w)^2 / sum(w^2); prev <- sum(w * y) / sum(w)
    return(list(constant = (neff * prev + .5) / (neff + 1), spec = spec,
                status = "sparse-class smoothed intercept fallback"))
  }
  e <- spec$engine
  model <- switch(e,
    glm = do.call(stats::glm.fit, c(list(x = cbind(1, x), y = y, weights = w, family = stats::quasibinomial()), a)),
    glmnet = {
      .need("glmnet"); a <- utils::modifyList(list(alpha = .5, lambda = .01), a)
      if (length(a$lambda) != 1L) stop("Each learner candidate must have one lambda.")
      do.call(glmnet::glmnet, c(list(x = x, y = y, weights = w, family = "binomial"), a))
    },
    ranger = {
      .need("ranger"); a <- utils::modifyList(list(num.trees = 500L, num.threads = 1L), a)
      do.call(ranger::ranger, c(list(x = as.data.frame(x), y = factor(y, levels = 0:1), case.weights = w, probability = TRUE), a))
    },
    xgboost = {
      .need("xgboost"); nr <- if (is.null(a$nrounds)) 100L else a$nrounds; a$nrounds <- NULL
      a <- utils::modifyList(list(objective = "binary:logistic", max_depth = 3L, eta = .05, nthread = 1L), a)
      xgboost::xgb.train(params = a, data = xgboost::xgb.DMatrix(x, label = y, weight = w), nrounds = nr, verbose = 0)
    },
    nnet = {
      .need("nnet"); a <- utils::modifyList(list(size = 3L, decay = .01, trace = FALSE, maxit = 200L), a)
      do.call(nnet::nnet, c(list(x = x, y = y, weights = w, entropy = TRUE), a))
    },
    custom = spec$fit(x, y, w))
  status <- "ok"
  if (e == "glmnet" && (!length(model$lambda) || (!is.null(model$jerr) && model$jerr != 0)))
    stop("glmnet did not return a converged candidate.")
  if (e == "nnet" && !is.null(model$convergence) && model$convergence != 0)
    stop("nnet did not converge; increase maxit or change the candidate.")
  if (e == "glm") {
    if (!model$converged) stop("GLM did not converge; use a regularized learner.")
    if (anyNA(model$coefficients)) {
      model$coefficients[is.na(model$coefficients)] <- 0
      status <- "rank-deficient columns omitted"
    }
  }
  # Portable bytes avoid XGBoost ALTREP deserialization before its namespace
  # has been loaded in a new session. Optional engines stay lazy at attachment.
  if (e == "xgboost") model <- xgboost::xgb.save.raw(model)
  list(model = model, spec = spec, status = status)
}
.predict_engine <- function(fit, x) {
  if (!is.null(fit$constant)) return(rep(.clip(fit$constant), nrow(x)))
  if (fit$spec$engine %in% c("glmnet", "ranger", "xgboost", "nnet")) .need(fit$spec$engine)
  p <- switch(fit$spec$engine,
    glm = stats::plogis(as.numeric(cbind(1, x) %*% fit$model$coefficients)),
    glmnet = as.numeric(stats::predict(fit$model, newx = x, type = "response", s = fit$model$lambda[1])),
    ranger = as.numeric(stats::predict(fit$model, data = as.data.frame(x))$predictions[, "1"]),
    xgboost = as.numeric(stats::predict(xgboost::xgb.load.raw(fit$model), x)),
    nnet = as.numeric(stats::predict(fit$model, x, type = "raw")),
    custom = as.numeric(fit$spec$predict(fit$model, x)))
  .prob(p, nrow(x))
}
.base_fit <- function(d, outcome, predictors, spec, recipe_args) {
  rec <- do.call(svy_recipe, c(list(data = d$data[predictors], w = d$w), recipe_args))
  fit <- .fit_engine(bake_svy_recipe(rec, d$data), .binary(d$data[[outcome]]), d$w, spec)
  list(recipe = rec, fit = fit)
}
.base_predict <- function(f, newdata) .predict_engine(f$fit, bake_svy_recipe(f$recipe, newdata))
.inner <- function(d, outcome, predictors, learners, folds, recipe_args) {
  y <- .binary(d$data[[outcome]]); mat <- matrix(NA_real_, nrow(d$data), length(learners))
  statuses <- rep("ok", length(learners))
  for (k in sort(unique(folds))) {
    tr <- .subset_design(d, which(folds != k)); te <- which(folds == k)
    for (j in seq_along(learners)) {
      tryCatch({
        f <- .base_fit(tr, outcome, predictors, learners[[j]], recipe_args)
        mat[te, j] <- .base_predict(f, d$data[te, , drop = FALSE])
        if (f$fit$status != "ok") statuses[j] <- if (statuses[j] == "ok") f$fit$status else
          paste(statuses[j], f$fit$status, sep = "; ")
      }, error = function(e) statuses[j] <<- conditionMessage(e))
    }
  }
  loss <- vapply(seq_along(learners), function(j)
    if (all(is.finite(mat[, j]))) svy_metrics(y, mat[, j], d$w)["log_loss"] else Inf, numeric(1))
  if (!any(is.finite(loss))) stop("All learner candidates failed inner validation: ", paste(unique(statuses), collapse = "; "))
  list(predictions = mat, loss = loss, status = statuses)
}
.stack_weights <- function(y, p, w) {
  if (ncol(p) == 1L) return(1)
  soft <- function(a) {z <- exp(c(a, 0) - max(c(a, 0))); z / sum(z)}
  opt <- stats::optim(rep(0, ncol(p) - 1L), function(a)
    svy_metrics(y, as.numeric(p %*% soft(a)), w)["log_loss"], method = "BFGS")
  if (opt$convergence != 0) stop("Stacking optimizer did not converge.")
  soft(opt$par)
}
.training_pipeline <- function(d, outcome, predictors, learners, inner_v, calibration, stack, recipe_args, seed) {
  # Isolate stochastic learners as well as fold allocation. An earlier outer
  # training set must not change this pipeline by consuming a different amount
  # of random state (for example after an outcome-dependent sparse fallback).
  .seed(seed, {
  nf <- min(inner_v, length(unique(.fold_cluster(d))))
  if (nf < 2L) stop("Too few training clusters for inner cross-fitting.")
  folds <- svy_split(d, nf, seed = seed)[, 1L]
  inn <- .inner(d, outcome, predictors, learners, folds, recipe_args)
  mix <- if (stack) {
    z <- numeric(length(learners)); ok <- is.finite(inn$loss)
    z[ok] <- .stack_weights(d$data[[outcome]], inn$predictions[, ok, drop = FALSE], d$w); z
  } else {
    z <- numeric(length(learners)); z[which.min(inn$loss)] <- 1; z
  }
  p <- as.numeric(inn$predictions[, mix > 0, drop = FALSE] %*% mix[mix > 0])
  cal <- svy_calibrate(d$data[[outcome]], p, d$w, method = calibration)
  fits <- lapply(which(mix > 0), function(j) .base_fit(d, outcome, predictors, learners[[j]], recipe_args))
  list(fits = fits, mixture = mix[mix > 0], candidate_weights = mix,
       calibrator = cal, inner_loss = inn$loss, inner_status = inn$status, inner_folds = folds)
  })
}
.pipeline_predict <- function(f, data, calibrated = TRUE) {
  mat <- vapply(f$fits, .base_predict, numeric(nrow(data)), newdata = data)
  mat <- matrix(mat, nrow(data), length(f$fits))
  p <- as.numeric(mat %*% f$mixture)
  if (calibrated) stats::predict(f$calibrator, p) else p
}

#' Nested whole-cluster learning and probability calibration
#'
#' Each outer validation fold is excluded from preprocessing, inner model
#' selection, stacking and calibration. Inner out-of-fold predictions train
#' the calibrator, and base models are refitted on the outer training sample.
#' All repeats are retained separately. Averaged predictions are explicitly
#' labeled as consensus; their performance differs from mean repeat performance.
#' @param design A svy_learn_design.
#' @param outcome Name of the outcome column; see family.
#' @param predictors Names of predictor columns; excludes design and outcome columns.
#' @param learners List of svy_learner specifications, treated as a tuning grid.
#' @param v,repeats,inner_v Outer folds, repeats and inner folds.
#' @param calibration Binary: platt, isotonic, spline or none. Other families:
#'   auto (default when omitted) or none. Auto uses affine mean calibration
#'   for gaussian, a mean ratio for poisson, and temperature scaling for multinomial.
#' @param stack Fit a nonnegative convex ensemble of candidates using inner OOF scores.
#' @param recipe_args Named arguments passed to svy_recipe.
#' @param seed Random seed; caller RNG state is restored.
#' @param family binomial (default), gaussian, poisson or multinomial.
#'   Gaussian accepts finite numeric outcomes; poisson accepts nonnegative
#'   integer counts. Multinomial requires a factor whose declared levels define
#'   the class schema before resampling. All outcomes must be observed.
#' @return A svy_learn_fit with OOF predictions, repeat metrics, fold assignments,
#'   selected candidate weights, final fitted pipeline and diagnostics.
#'   Nonbinary fits additionally inherit svy_response_fit and record family.
#'   Scalar predictions have dimensions observations by repeats; multinomial
#'   predictions have dimensions observations by classes by repeats. Consensus
#'   predictions average the repeats. Nonbinary calibration can be disabled by
#'   calibration = "none"; calibrated = FALSE in predict returns raw predictions.
#' @param attrition Optional unfitted svy_attrition specification. Supply a full
#'   panel roster and baseline weights. Retention models are fitted separately in
#'   every inner and outer training partition; outcome fitting uses observed,
#'   eligible rows. Results retain repeat-specific cross-fitted weights and the
#'   baseline roster for complete-workflow replication.
#' @export
#' @examples
#' d <- svy_example(200)
#' design <- svy_design(d, "weight", "psu", "stratum")
#' fit <- svy_fit(design, "y", c("x1", "x2"), v = 3, inner_v = 2)
#' summary(fit)
svy_fit <- function(design, outcome, predictors, learners = list(svy_learner()),
                    v = 5L, repeats = 1L, inner_v = 3L, calibration = "platt",
                    stack = FALSE, recipe_args = list(), seed = 1L,
                    family = "binomial", attrition = NULL) {
  family <- match.arg(family, c("binomial", "gaussian", "poisson", "multinomial"))
  if (!is.null(attrition)) return(.attrition_fit(design, outcome, predictors,
    learners, v, repeats, inner_v,
    if (family != "binomial" && missing(calibration)) "auto" else calibration,
    stack, recipe_args, seed, family, attrition, match.call()))
  if (family != "binomial") return(.response_fit(design, outcome, predictors,
    learners, v, repeats, inner_v, if (missing(calibration)) "auto" else calibration,
    stack, recipe_args, seed, family, match.call()))
  y <- .binary(design$data[[outcome]])
  if (!length(predictors) || !all(predictors %in% names(design$data)) ||
      any(predictors %in% c(outcome, design$weights, design$ids, design$strata,
                            design$fpc, design$subject, design$household, design$fold_groups)))
    stop("Predictors must exist and exclude outcome, weights and sampling identifiers.")
  if (!length(learners) || !all(vapply(learners, inherits, logical(1), "svy_learner")))
    stop("learners must be a nonempty list of svy_learner objects.")
  if (inner_v < 2L || inner_v != as.integer(inner_v)) stop("inner_v must be an integer >= 2.")
  .seed(seed, {
    folds <- svy_split(design, v, repeats, seed)
    raw <- calibrated <- matrix(NA_real_, nrow(design$data), repeats)
    audit <- list(); outer_fits <- list()
    for (r in seq_len(repeats)) for (k in seq_len(v)) {
      tr <- which(folds[, r] != k); te <- which(folds[, r] == k)
      f <- .training_pipeline(.subset_design(design, tr), outcome, predictors,
        learners, inner_v, calibration, stack, recipe_args, seed + r * 1000L + k)
      raw[te, r] <- .pipeline_predict(f, design$data[te, , drop = FALSE], FALSE)
      calibrated[te, r] <- stats::predict(f$calibrator, raw[te, r])
      audit[[length(audit) + 1L]] <- data.frame(repeat_id = r, fold = k,
        train_n = length(tr), test_n = length(te), train_events = sum(y[tr]),
        test_events = sum(y[te]), candidate = seq_along(learners),
        inner_log_loss = f$inner_loss, mixture_weight = f$candidate_weights,
        inner_status = f$inner_status,
        calibration_status = f$calibrator$status,
        learner_status = paste(vapply(f$fits, function(z) z$fit$status, character(1)), collapse = "; "))
      outer_fits[[length(outer_fits) + 1L]] <- f
    }
    metrics <- do.call(rbind, lapply(seq_len(repeats), function(r)
      svy_metrics(y, calibrated[, r], design$w)))
    final <- .training_pipeline(design, outcome, predictors, learners, inner_v,
      calibration, stack, recipe_args, seed + 900000L)
    structure(list(final = final, outer_fits = outer_fits, folds = folds,
      predictions_raw = raw, predictions = calibrated, consensus = rowMeans(calibrated),
      metrics = metrics, audit = do.call(rbind, audit), outcome = outcome,
      predictors = predictors, design = design, diagnostics = svy_diagnostics(design, y),
      configuration = list(outcome = outcome, predictors = predictors, learners = learners,
        v = v, repeats = repeats, inner_v = inner_v, calibration = calibration,
        stack = stack, recipe_args = recipe_args, seed = seed),
      call = match.call(), seed = seed), class = "svy_learn_fit")
  })
}

#' @export
predict.svy_learn_fit <- function(object, newdata, calibrated = TRUE, ...) {
  .pipeline_predict(object$final, newdata, calibrated)
}
#' @export
summary.svy_learn_fit <- function(object, ...) {
  list(design = object$diagnostics,
    mean_repeat_metrics = colMeans(object$metrics, na.rm = TRUE),
    repeat_sd = apply(object$metrics, 2L, stats::sd),
    consensus_metrics = svy_metrics(object$design$data[[object$outcome]], object$consensus, object$design$w),
    final_candidate_weights = object$final$candidate_weights,
    calibration_status = unique(object$audit$calibration_status))
}
#' @export
print.svy_learn_fit <- function(x, ...) {
  cat("Survey-aware binary learning\n")
  cat(x$diagnostics$n, "observations;", x$diagnostics$psus, "PSUs;",
      ncol(x$predictions), "CV repeats\n")
  print(colMeans(x$metrics, na.rm = TRUE)); invisible(x)
}
#' @export
plot.svy_learn_fit <- function(x, ...) {
  d <- svy_calibration_bins(x$design$data[[x$outcome]], x$consensus, x$design$w)
  graphics::plot(d$predicted, d$observed, xlim = c(0, 1), ylim = c(0, 1),
    xlab = "Weighted mean predicted probability", ylab = "Weighted observed proportion", ...)
  graphics::abline(0, 1, lty = 2, col = "gray50"); invisible(d)
}
