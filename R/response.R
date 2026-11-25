# Outcome-specific operations share the design, frozen recipe, cluster splitting
# and nested evaluation contract with the binary workflow.
.response_y <- function(y, family) {
  if (family == "multinomial") {
    if (!is.factor(y) || nlevels(y) < 2L || !length(y) || anyNA(y))
      stop("Multinomial outcomes require a nonmissing factor with declared class levels.")
  } else {
    if (!is.numeric(y) || !length(y) || any(!is.finite(y)))
      stop("Outcomes must be a nonempty finite numeric vector.")
    if (family == "poisson" && any(y < 0 | y != floor(y)))
      stop("Poisson outcomes must be nonnegative integer counts.")
  }
  y
}
.response_p <- function(p, n, family, levels = NULL) {
  if (family == "multinomial") {
    if (!is.matrix(p) || !is.numeric(p) || nrow(p) != n || ncol(p) != length(levels) ||
        is.null(colnames(p)) || anyDuplicated(colnames(p)) || !setequal(colnames(p), levels) ||
        any(!is.finite(p)) || any(p < 0 | p > 1) || any(abs(rowSums(p) - 1) > 1e-6))
      stop("Class probabilities must be a finite row-normalized matrix with columns matching declared levels.")
    p <- p[, levels, drop = FALSE]
    return(p / rowSums(p))
  }
  if (!is.numeric(p) || length(p) != n || any(!is.finite(p)) ||
      (family == "poisson" && any(p < 0))) stop("Invalid response predictions.")
  as.numeric(p)
}
.response_losses <- function(y, p, family) {
  if (family == "multinomial") {
    hit <- p[cbind(seq_along(y), as.integer(y))]
    return(cbind(log_loss = -log(pmax(hit, 1e-15)),
      brier = rowSums(p^2) - 2 * hit + 1,
      accuracy = as.numeric(max.col(p, ties.method = "first") == as.integer(y))))
  }
  out <- cbind(mse = (y - p)^2, mae = abs(y - p), bias = p - y)
  if (family == "poisson") {
    mu <- pmax(p, 1e-15); term <- numeric(length(y)); positive <- y > 0
    term[positive] <- y[positive] * (log(y[positive]) - log(mu[positive]))
    out <- cbind(poisson_deviance = pmax(0, 2 * (term - y + mu)), out)
  }
  out
}
.response_metrics <- function(y, p, w, family) {
  y <- .response_y(y, family); p <- .response_p(p, length(y), family, levels(y))
  w <- .positive(w, length(y), TRUE); w <- w / sum(w)
  out <- colSums(.response_losses(y, p, family) * w)
  if (family != "multinomial") out <- c(out, rmse = sqrt(out[["mse"]]),
    mean_outcome = sum(w * y), mean_prediction = sum(w * p))
  out
}
.softmax_rows <- function(z) {
  z <- z - apply(z, 1L, max)
  e <- exp(z); e / rowSums(e)
}
.response_engine <- function(x, y, w, spec, family) {
  w <- w / mean(w); e <- spec$engine; a <- spec$parameters
  lev <- levels(y); status <- "ok"
  out <- list(spec = spec, family = family, levels = lev, status = status)
  if (family == "multinomial" && min(tabulate(as.integer(y), length(lev))) < 2L) {
    neff <- sum(w)^2 / sum(w^2)
    freq <- vapply(lev, function(k) sum(w[y == k]) / sum(w), numeric(1))
    out$constant <- (neff * freq + .5) / (neff + length(lev) * .5)
    out$status <- "sparse-class smoothed intercept fallback"; return(out)
  }
  if (family != "multinomial" && length(unique(y)) == 1L) {
    out$constant <- y[1L]; out$status <- "constant response fallback"; return(out)
  }
  model <- switch(e,
    glm = {
      if (family == "multinomial") {
        if (length(setdiff(names(a), c("lambda", "maxit")))) stop("Multinomial glm accepts lambda and maxit only.")
        lambda <- if (is.null(a$lambda)) 1e-6 else a$lambda
        if (length(lambda) != 1L || !is.finite(lambda) || lambda < 0) stop("Invalid multinomial ridge penalty.")
        xx <- cbind(1, x); k <- length(lev); q <- ncol(xx)
        target <- matrix(0, length(y), k); target[cbind(seq_along(y), as.integer(y))] <- 1
        ww <- w / sum(w)
        value <- function(b, gradient = FALSE) {
          bb <- matrix(b, q, k - 1L); eta <- cbind(xx %*% bb, 0)
          m <- apply(eta, 1L, max); logden <- m + log(rowSums(exp(eta - m)))
          if (!gradient) return(sum(ww * (logden - eta[cbind(seq_along(y), as.integer(y))])) +
            lambda * sum(bb[-1L, , drop = FALSE]^2) / 2)
          grad <- crossprod(xx, (.softmax_rows(eta) - target)[, seq_len(k - 1L), drop = FALSE] * ww)
          grad[-1L, ] <- grad[-1L, , drop = FALSE] + lambda * bb[-1L, , drop = FALSE]
          as.numeric(grad)
        }
        opt <- stats::optim(rep(0, q * (k - 1L)), value,
          function(b) value(b, TRUE), method = "BFGS",
          control = list(maxit = if (is.null(a$maxit)) 5000L else a$maxit, reltol = 1e-10))
        if (opt$convergence != 0) stop("Multinomial GLM did not converge.")
        list(coefficients = matrix(opt$par, q, k - 1L))
      } else {
        f <- do.call(stats::glm.fit, c(list(x = cbind(1, x), y = y, weights = w,
          family = if (family == "gaussian") stats::gaussian() else stats::poisson()), a))
        if (!f$converged) stop("GLM did not converge.")
        if (anyNA(f$coefficients)) {
          f$coefficients[is.na(f$coefficients)] <- 0
          out$status <- "rank-deficient columns omitted"
        }
        f
      }
    },
    glmnet = {
      .need("glmnet"); a <- utils::modifyList(list(alpha = .5, lambda = .01), a)
      if (length(a$lambda) != 1L) stop("Each learner candidate must have one lambda.")
      f <- do.call(glmnet::glmnet, c(list(x = x, y = y, weights = w, family = family), a))
      if (!length(f$lambda) || (!is.null(f$jerr) && f$jerr != 0)) stop("glmnet did not return a converged candidate.")
      f
    },
    ranger = {
      .need("ranger"); a <- utils::modifyList(list(num.trees = 500L, num.threads = 1L), a)
      do.call(ranger::ranger, c(list(x = as.data.frame(x), y = y, case.weights = w,
        probability = family == "multinomial"), a))
    },
    xgboost = {
      .need("xgboost"); nr <- if (is.null(a$nrounds)) 100L else a$nrounds; a$nrounds <- NULL
      defaults <- list(objective = switch(family, gaussian = "reg:squarederror",
        poisson = "count:poisson", multinomial = "multi:softprob"),
        max_depth = 3L, eta = .05, nthread = 1L)
      if (family == "multinomial") defaults$num_class <- length(lev)
      a <- utils::modifyList(defaults, a)
      yy <- if (family == "multinomial") as.integer(y) - 1L else y
      xgboost::xgb.save.raw(xgboost::xgb.train(params = a,
        data = xgboost::xgb.DMatrix(x, label = yy, weight = w), nrounds = nr, verbose = 0))
    },
    nnet = {
      .need("nnet")
      if (family == "poisson") stop("The nnet adapter supports gaussian and multinomial responses; use glm, glmnet, ranger or xgboost for counts.")
      a <- utils::modifyList(list(size = 3L, decay = .01, trace = FALSE, maxit = 500L), a)
      yy <- if (family == "multinomial") {
        z <- matrix(0, length(y), length(lev)); z[cbind(seq_along(y), as.integer(y))] <- 1; z
      } else y
      f <- do.call(nnet::nnet, c(list(x = x, y = yy, weights = w,
        linout = family == "gaussian", softmax = family == "multinomial"), a))
      if (f$convergence != 0) stop("nnet did not converge.")
      f
    },
    custom = spec$fit(x, y, w))
  out$model <- model; out
}
.response_predict_engine <- function(f, x) {
  n <- nrow(x); k <- length(f$levels); family <- f$family
  if (!is.null(f$constant)) {
    if (family == "multinomial") return(matrix(rep(f$constant, each = n), n, k,
      dimnames = list(NULL, f$levels)))
    return(rep(f$constant, n))
  }
  e <- f$spec$engine
  if (e %in% c("glmnet", "ranger", "xgboost", "nnet")) .need(e)
  p <- switch(e,
    glm = {
      eta <- cbind(1, x) %*% f$model$coefficients
      if (family == "multinomial") .softmax_rows(cbind(eta, 0)) else
        if (family == "poisson") exp(as.numeric(eta)) else as.numeric(eta)
    },
    glmnet = stats::predict(f$model, newx = x, type = "response", s = f$model$lambda[1L]),
    ranger = stats::predict(f$model, data = as.data.frame(x))$predictions,
    xgboost = stats::predict(xgboost::xgb.load.raw(f$model), x),
    nnet = stats::predict(f$model, x, type = "raw"),
    custom = f$spec$predict(f$model, x))
  if (family == "multinomial" && e != "custom") {
    # glmnet supplies an n x class x lambda array; XGBoost may supply row-major values.
    labels <- if (e == "glmnet") dimnames(p)[[2L]] else colnames(p)
    if (e == "xgboost" && is.null(dim(p))) p <- matrix(p, n, k, byrow = TRUE) else p <- matrix(p, n, k)
    colnames(p) <- if (is.null(labels)) f$levels else labels
  }
  .response_p(p, n, family, f$levels)
}
.response_base <- function(d, outcome, predictors, spec, recipe_args, family) {
  rec <- do.call(svy_recipe, c(list(data = d$data[predictors], w = d$w), recipe_args))
  list(recipe = rec, fit = .response_engine(bake_svy_recipe(rec, d$data),
    d$data[[outcome]], d$w, spec, family))
}
.response_base_predict <- function(f, data) .response_predict_engine(f$fit, bake_svy_recipe(f$recipe, data))
.response_calibrator <- function(y, p, w, family, method) {
  out <- list(family = family, method = method, status = "ok")
  if (method == "none") return(out)
  if (family == "gaussian") {
    m <- stats::lm.wfit(cbind(1, p), y, w / mean(w))
    b <- m$coefficients
    if (anyNA(b)) { b <- c(sum(w * y) / sum(w), 0); out$status <- "constant calibration fallback" }
    out$coefficients <- b
  } else if (family == "poisson") {
    den <- sum(w * p)
    if (den <= 1e-15 * sum(w)) {
      out$constant <- sum(w * y) / sum(w); out$status <- "constant calibration fallback"
    } else out$ratio <- sum(w * y) / den
  } else {
    # A bounded positive temperature changes sharpness and preserves class ranking.
    loss <- function(logt) .response_metrics(y,
      .response_temperature(p, exp(logt)), w, family)[["log_loss"]]
    opt <- stats::optimize(loss, interval = c(-4, 4), tol = 1e-8)
    out$temperature <- exp(opt$minimum)
    if (abs(opt$minimum) > 3.99) out$status <- "temperature at search boundary"
  }
  out
}
.response_temperature <- function(p, temperature) {
  out <- .softmax_rows(log(pmax(p, 1e-15)) / temperature)
  colnames(out) <- colnames(p); out
}
.response_calibrate_predict <- function(cal, p) {
  if (cal$method == "none") return(p)
  if (cal$family == "gaussian") return(as.numeric(cal$coefficients[1] + cal$coefficients[2] * p))
  if (cal$family == "poisson") return(if (!is.null(cal$constant)) rep(cal$constant, length(p)) else p * cal$ratio)
  .response_temperature(p, cal$temperature)
}
.response_mix <- function(predictions, mix) {
  active <- which(mix > 0)
  Reduce(`+`, lapply(active, function(j) predictions[[j]] * mix[j]))
}
.response_pipeline <- function(d, outcome, predictors, learners, inner_v, calibration,
                               stack, recipe_args, seed, family) .seed(seed, {
  nf <- min(inner_v, length(unique(.fold_cluster(d))))
  if (nf < 2L) stop("Too few training clusters for inner cross-fitting.")
  folds <- svy_split(d, nf, seed = seed)[, 1L]
  y <- d$data[[outcome]]; n <- length(y); lev <- levels(y)
  blank <- if (family == "multinomial") matrix(NA_real_, n, length(lev), dimnames = list(NULL, lev)) else rep(NA_real_, n)
  preds <- rep(list(blank), length(learners)); status <- rep("ok", length(learners))
  failed <- rep(FALSE, length(learners))
  for (k in seq_len(nf)) {
    tr <- .subset_design(d, which(folds != k)); te <- which(folds == k)
    for (j in seq_along(learners)) tryCatch({
      f <- .response_base(tr, outcome, predictors, learners[[j]], recipe_args, family)
      p <- .response_base_predict(f, d$data[te, , drop = FALSE])
      if (family == "multinomial") preds[[j]][te, ] <- p else preds[[j]][te] <- p
      if (f$fit$status != "ok") status[j] <- paste(status[j], f$fit$status, sep = "; ")
    }, error = function(e) { failed[j] <<- TRUE; status[j] <<- conditionMessage(e) })
  }
  loss <- vapply(seq_along(learners), function(j) if (failed[j]) Inf else
    .response_metrics(y, preds[[j]], d$w, family)[[1L]], numeric(1))
  good <- which(is.finite(loss))
  if (!length(good)) stop("All learner candidates failed inner validation: ", paste(unique(status), collapse = "; "))
  mix <- numeric(length(learners))
  if (stack && length(good) > 1L) {
    soft <- function(a) {z <- exp(c(a, 0) - max(c(a, 0))); z / sum(z)}
    opt <- stats::optim(rep(0, length(good) - 1L), function(a)
      .response_metrics(y, .response_mix(preds[good], soft(a)), d$w, family)[[1L]], method = "BFGS")
    if (opt$convergence != 0) stop("Stacking optimizer did not converge.")
    mix[good] <- soft(opt$par)
  } else mix[which.min(loss)] <- 1
  cal <- .response_calibrator(y, .response_mix(preds, mix), d$w, family, calibration)
  fits <- lapply(which(mix > 0), function(j) .response_base(d, outcome, predictors,
    learners[[j]], recipe_args, family))
  list(fits = fits, mixture = mix[mix > 0], candidate_weights = mix, calibrator = cal,
    inner_loss = loss, inner_status = status, inner_folds = folds)
})
.response_pipeline_predict <- function(f, data, calibrated = TRUE) {
  p <- .response_mix(lapply(f$fits, .response_base_predict, data = data), f$mixture)
  if (calibrated) .response_calibrate_predict(f$calibrator, p) else p
}
.response_fit <- function(design, outcome, predictors, learners, v, repeats, inner_v,
                         calibration, stack, recipe_args, seed, family, call) {
  y <- .response_y(design$data[[outcome]], family)
  if (!length(predictors) || !all(predictors %in% names(design$data)) ||
      any(predictors %in% c(outcome, design$weights, design$ids, design$strata,
        design$fpc, design$subject, design$household, design$fold_groups))) stop("Predictors must exist and exclude outcome, weights and sampling identifiers.")
  if (!length(learners) || !all(vapply(learners, inherits, logical(1), "svy_learner")))
    stop("learners must be a nonempty list of svy_learner objects.")
  if (length(inner_v) != 1L || !is.finite(inner_v) || inner_v < 2L || inner_v != as.integer(inner_v))
    stop("inner_v must be an integer >= 2.")
  calibration <- match.arg(calibration, c("auto", "none"))
  folds <- svy_split(design, v, repeats, seed)
  n <- length(y); lev <- levels(y); categorical <- family == "multinomial"
  raw <- pred <- if (categorical) array(NA_real_, c(n, length(lev), repeats),
    dimnames = list(NULL, lev, NULL)) else matrix(NA_real_, n, repeats)
  audit <- outer <- list()
  for (r in seq_len(repeats)) for (k in seq_len(v)) {
    tr <- which(folds[, r] != k); te <- which(folds[, r] == k)
    f <- .response_pipeline(.subset_design(design, tr), outcome, predictors, learners,
      inner_v, calibration, stack, recipe_args, seed + r * 1000L + k, family)
    p <- .response_pipeline_predict(f, design$data[te, , drop = FALSE], FALSE)
    pp <- .response_calibrate_predict(f$calibrator, p)
    if (categorical) { raw[te, , r] <- p; pred[te, , r] <- pp } else { raw[te, r] <- p; pred[te, r] <- pp }
    audit[[length(audit) + 1L]] <- data.frame(repeat_id = r, fold = k,
      train_n = length(tr), test_n = length(te), candidate = seq_along(learners),
      inner_loss = f$inner_loss, mixture_weight = f$candidate_weights, inner_status = f$inner_status,
      calibration_status = f$calibrator$status,
      learner_status = paste(vapply(f$fits, function(z) z$fit$status, character(1)), collapse = "; "))
    outer[[length(outer) + 1L]] <- f
  }
  metrics <- do.call(rbind, lapply(seq_len(repeats), function(r) .response_metrics(y,
    if (categorical) matrix(pred[, , r], n, length(lev), dimnames = list(NULL, lev)) else pred[, r], design$w, family)))
  consensus <- if (categorical) apply(pred, c(1L, 2L), mean) else rowMeans(pred)
  final <- .response_pipeline(design, outcome, predictors, learners, inner_v,
    calibration, stack, recipe_args, seed + 900000L, family)
  structure(list(final = final, outer_fits = outer, folds = folds, predictions_raw = raw,
    predictions = pred, consensus = consensus, metrics = metrics, audit = do.call(rbind, audit),
    outcome = outcome, predictors = predictors, family = family, levels = lev, design = design,
    diagnostics = svy_diagnostics(design), configuration = list(outcome = outcome, predictors = predictors,
      learners = learners, v = v, repeats = repeats, inner_v = inner_v, calibration = calibration,
      stack = stack, recipe_args = recipe_args, seed = seed, family = family), call = call, seed = seed),
    class = c("svy_response_fit", "svy_learn_fit"))
}

#' @export
predict.svy_response_fit <- function(object, newdata, calibrated = TRUE, ...) {
  .response_pipeline_predict(object$final, newdata, calibrated)
}
#' @export
summary.svy_response_fit <- function(object, ...) {
  list(family = object$family, design = object$diagnostics,
    mean_repeat_metrics = colMeans(object$metrics), repeat_sd = apply(object$metrics, 2L, stats::sd),
    consensus_metrics = svy_metrics(object$design$data[[object$outcome]], object$consensus,
      object$design$w, object$family), final_candidate_weights = object$final$candidate_weights,
    calibration_status = unique(object$audit$calibration_status))
}
#' @export
print.svy_response_fit <- function(x, ...) {
  cat("Survey-aware", x$family, "learning\n")
  cat(x$diagnostics$n, "observations;", x$diagnostics$psus, "PSUs;", ncol(x$folds), "CV repeats\n")
  print(colMeans(x$metrics)); invisible(x)
}
#' @export
plot.svy_response_fit <- function(x, ...) {
  y <- x$design$data[[x$outcome]]; p <- x$consensus; w <- x$design$w
  if (x$family == "multinomial") {
    observed <- vapply(x$levels, function(k) sum(w * (y == k)) / sum(w), numeric(1))
    predicted <- colSums(p * w) / sum(w)
    graphics::plot(predicted, observed, xlim = c(0, 1), ylim = c(0, 1),
      xlab = "Weighted mean class probability", ylab = "Weighted class proportion", ...)
    graphics::text(predicted, observed, labels = x$levels, pos = 3)
  } else graphics::plot(p, y, xlab = "Out-of-fold predicted response", ylab = "Observed response", ...)
  graphics::abline(0, 1, lty = 2, col = "gray50"); invisible(x)
}
.response_infer <- function(design, y, p, domain, level, replicate_design, family) {
  .need("survey"); y <- .response_y(y, family)
  p <- .response_p(p, length(y), family, levels(y))
  if (length(y) != nrow(design$data)) stop("Row count mismatch.")
  if (!is.logical(domain) || length(domain) != length(y) || anyNA(domain) || !any(domain)) stop("Invalid domain.")
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) stop("Invalid confidence level.")
  opt <- options(survey.lonely.psu = design$lonely.psu); on.exit(options(opt), add = TRUE)
  des <- if (is.null(replicate_design)) as_survey_design(design) else replicate_design
  if (nrow(des$variables) != length(y)) stop("Replicate-design row count mismatch.")
  losses <- .response_losses(y, p, family)
  vars <- paste0(".svyl_response_", colnames(losses))
  for (j in seq_along(vars)) des$variables[[vars[j]]] <- losses[, j]
  des <- .survey_domain(des, domain); df <- survey::degf(des)
  if (df <= 0) stop("Positive survey degrees of freedom are required for intervals.")
  est <- survey::svymean(stats::reformulate(vars), des, deff = "replace")
  val <- as.numeric(stats::coef(est)); se <- sqrt(diag(stats::vcov(est)))
  q <- stats::qt((1 + level) / 2, df)
  scores <- data.frame(metric = colnames(losses), estimate = val, se = se,
    lower = val - q * se, upper = val + q * se, deff = as.numeric(survey::deff(est)), df = df, row.names = NULL)
  list(scores = scores, calibration = NULL, ranks = NULL,
    metrics = svy_metrics(y[domain], if (is.matrix(p)) p[domain, , drop = FALSE] else p[domain],
      design$w[domain], family), target = "Conditional mean-loss uncertainty for the supplied predictions")
}
