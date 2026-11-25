#' Weighted split-conformal thresholds and prediction sets
#'
#' Calibration predictions must come from a model fitted independently of the
#' calibration outcomes. Pooled cross-validation residuals are not split conformal.
#' The test point contributes an atom at positive infinity. Exchangeable mode
#' requires equal weights and exchangeable calibration/test units. Covariate-shift
#' mode requires the appropriate density ratios and invariant conditional outcome
#' distribution. Survey mode is an empirical weighted construction with no automatic
#' finite-sample coverage guarantee for dependent or informatively sampled units.
#' @param y,p Calibration outcomes and predictions.
#' @param newdata Test predictions, not a predictor data frame.
#' @param w Calibration weights or density ratios, according to weight_type.
#' @param test_weights Test weights on the same scale as w.
#' @param alpha Miscoverage level in (0,1).
#' @param family gaussian, binomial or multinomial. Gaussian uses absolute residuals;
#'   classification uses one minus the probability of the observed class.
#' @param weight_type exchangeable, covariate_shift or survey.
#' @return A svy_conformal object with thresholds and intervals or membership matrix.
#' @export
svy_conformal <- function(y, p, newdata, w = rep(1, length(y)),
                          test_weights = NULL, alpha = .1, family = "gaussian",
                          weight_type = "exchangeable") {
  family <- match.arg(family, c("gaussian", "binomial", "multinomial"))
  weight_type <- match.arg(weight_type, c("exchangeable", "covariate_shift", "survey"))
  .eval_level(alpha); w <- .positive(w, length(y), TRUE)
  if (family == "binomial") {
    y <- .binary(y); p <- .prob(p, length(y)); pp <- .prob(newdata, length(newdata))
    scores <- ifelse(y == 1, 1 - p, p)
    pred <- cbind(`0` = 1 - pp, `1` = pp)
  } else if (family == "multinomial") {
    y <- .response_y(y, family); p <- .response_p(p, length(y), family, levels(y))
    pred <- .response_p(newdata, nrow(newdata), family, levels(y))
    scores <- 1 - p[cbind(seq_along(y), as.integer(y))]
  } else {
    y <- .response_y(y, family); p <- .response_p(p, length(y), family)
    pred <- .response_p(newdata, length(newdata), family); scores <- abs(y - p)
  }
  nt <- if (is.matrix(pred)) nrow(pred) else length(pred)
  if (is.null(test_weights)) {
    if (weight_type != "exchangeable") stop("Supply test_weights on the calibration-weight scale.")
    test_weights <- rep(w[1L], nt)
  }
  test_weights <- .positive(test_weights, nt)
  if (weight_type == "exchangeable" &&
      (any(w != w[1L]) || any(test_weights != w[1L])))
    stop("Exchangeable mode requires identical calibration and test weights.")
  ok <- w > 0; o <- order(scores[ok]); sc <- scores[ok][o]; ww <- w[ok][o]
  # Normalize before summation to avoid overflow and preserve common-scale invariance.
  scale <- max(c(ww, test_weights)); ww <- ww / scale; tw <- test_weights / scale
  cw <- cumsum(ww); threshold_mass <- (1 - alpha) * (sum(ww) + tw)
  # findInterval(left.open=TRUE) gives the first cumulative mass >= target, including ties.
  k <- findInterval(threshold_mass, cw, left.open = TRUE) + 1L
  q <- rep(Inf, nt); finite <- k <= length(sc); q[finite] <- sc[k[finite]]
  sets <- if (family == "gaussian") cbind(lower = pred - q, upper = pred + q) else (1 - pred) <= q
  structure(list(sets = sets, threshold = q, alpha = alpha, family = family,
    weight_type = weight_type, calibration_n = length(y),
    assumption = switch(weight_type,
      exchangeable = "Independent fitted rule; exchangeable calibration and test units; marginal coverage",
      covariate_shift = "Independent fitted rule; correct density ratios; invariant conditional outcome law; marginal coverage",
      survey = "Empirical survey weighting; finite-sample coverage is not guaranteed")), class = "svy_conformal")
}

#' Evaluate prediction sets on a separate probability sample
#' @param design A svy_learn_design for test observations only.
#' @param y Test outcomes.
#' @param sets A svy_conformal object constructed without the test outcomes.
#' @param level Confidence level for conditional design-based summaries.
#' @param replicate_design Optional matching survey replicate design.
#' @return Empirical conf_coverage, conf_width and conditional intervals. Infinite
#'   interval width is retained with an undefined width standard error.
#' @export
svy_conformal_eval <- function(design, y, sets, level = .95, replicate_design = NULL) {
  .eval_level(level)
  if (!inherits(sets, "svy_conformal") || nrow(sets$sets) != length(y) ||
      length(y) != nrow(design$data)) stop("Prediction sets and test rows must match.")
  if (sets$family == "gaussian") {
    y <- .response_y(y, "gaussian")
    coverage <- y >= sets$sets[, 1L] & y <= sets$sets[, 2L]
    width <- sets$sets[, 2L] - sets$sets[, 1L]
  } else {
    lab <- if (sets$family == "binomial") as.character(.binary(y)) else as.character(y)
    if (anyNA(lab) || !all(lab %in% colnames(sets$sets))) stop("Unknown test outcome class.")
    coverage <- sets$sets[cbind(seq_along(y), match(lab, colnames(sets$sets)))]
    width <- rowSums(sets$sets)
  }
  opt <- options(survey.lonely.psu = design$lonely.psu); on.exit(options(opt), add = TRUE)
  des <- .eval_design(design, replicate_design)
  ci <- .eval_mean(des, cbind(conf_coverage = as.numeric(coverage)), level)
  if (all(is.finite(width))) ci <- rbind(ci, .eval_mean(des, cbind(conf_width = width), level)) else
    ci <- rbind(ci, data.frame(metric = "conf_width", estimate = Inf, se = NA_real_,
      ci_lower = NA_real_, ci_upper = NA_real_, df = survey::degf(des), method = "infinite width"))
  list(conf_coverage = sum(design$w * coverage) / sum(design$w),
    conf_width = sum(design$w * width) / sum(design$w), intervals = ci,
    target = "Empirical test-set coverage conditional on supplied sets", assumption = sets$assumption)
}
