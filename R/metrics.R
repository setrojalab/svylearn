#' Weighted prediction metrics
#'
#' AUC gives half credit for tied scores. Average precision integrates the
#' stepwise precision-recall curve at distinct score thresholds. Log loss clips
#' probabilities at 1e-15; the Brier score uses the original probabilities.
#' Gaussian scores include mean squared error, mean absolute error, signed
#' prediction bias and root mean squared error. Poisson scores additionally
#' include mean Poisson deviance, with mean predictions floored at 1e-15 for
#' scoring. Multinomial scores are log loss, the sum of squared class-probability
#' errors (without division by class count), and accuracy. Accuracy ties select
#' the first declared factor level. Probability row sums within 1e-6 of one
#' are normalized to remove roundoff before multiclass scoring.
#' @param y Outcomes: 0/1, numeric responses, counts or a factor, according to family.
#' @param p Predictions: numeric vector, or a matrix of class probabilities with
#'   named columns matching factor levels for multinomial outcomes.
#' @param w Nonnegative evaluation weights; defaults to equal weights.
#' @param family binomial, gaussian, poisson or multinomial.
#' @return Named numeric vector. The first metric is the family-specific
#'   selection loss: log_loss, mse, poisson_deviance or log_loss, respectively.
#'   Binary rank metrics are NA for a single outcome class.
#' @export
#' @examples
#' svy_metrics(c(0, 1, 0, 1), c(.1, .8, .2, .6))
svy_metrics <- function(y, p, w = rep(1, length(y)), family = "binomial") {
  family <- match.arg(family, c("binomial", "gaussian", "poisson", "multinomial"))
  if (family != "binomial") return(.response_metrics(y, p, w, family))
  y <- .binary(y); p <- .prob(p, length(y)); w <- .positive(w, length(y), zero = TRUE)
  w <- w / sum(w); prev <- sum(w * y); pc <- .clip(p)
  brier <- sum(w * (y - p)^2)
  auc <- ap <- NA_real_
  if (prev > 0 && prev < 1) {
    o <- order(p); ps <- p[o]; ys <- y[o]; ws <- w[o]
    g <- cumsum(c(TRUE, diff(ps) != 0))
    pos <- as.numeric(rowsum(ws * ys, g, reorder = FALSE))
    neg <- as.numeric(rowsum(ws * (1 - ys), g, reorder = FALSE))
    auc <- sum(pos * (cumsum(neg) - neg / 2)) / (prev * (1 - prev))
    pos <- rev(pos); neg <- rev(neg)
    denom <- cumsum(pos + neg)
    precision <- ifelse(denom > 0, cumsum(pos) / denom, 0)
    ap <- sum(pos / prev * precision)
  }
  c(log_loss = -sum(w * (y * log(pc) + (1 - y) * log1p(-pc))),
    brier = brier, scaled_brier = if (prev > 0 && prev < 1) 1 - brier / (prev * (1 - prev)) else NA_real_,
    auc = auc, average_precision = ap, prevalence = prev, mean_prediction = sum(w * p))
}

#' Weighted threshold and decision-curve summaries
#' @param y,p,w Outcomes, probabilities and nonnegative weights.
#' @param thresholds Decision thresholds strictly between zero and one.
#' @return A data frame including weighted confusion totals and net benefit.
#' @export
svy_decision <- function(y, p, w = rep(1, length(y)), thresholds = seq(.01, .3, .01)) {
  y <- .binary(y); p <- .prob(p, length(y)); w <- .positive(w, length(y), TRUE)
  if (!length(thresholds) || any(!is.finite(thresholds)) || any(thresholds <= 0 | thresholds >= 1))
    stop("thresholds must be strictly between zero and one.")
  ratio <- function(a, b) if (b > 0) a / b else NA_real_
  do.call(rbind, lapply(thresholds, function(t) {
    yes <- p >= t
    tp <- sum(w * (yes & y == 1)); fp <- sum(w * (yes & y == 0))
    tn <- sum(w * (!yes & y == 0)); fn <- sum(w * (!yes & y == 1))
    data.frame(threshold = t, tp = tp, fp = fp, tn = tn, fn = fn,
      sensitivity = ratio(tp, tp + fn), specificity = ratio(tn, tn + fp),
      ppv = ratio(tp, tp + fp), npv = ratio(tn, tn + fn),
      net_benefit = (tp - fp * t / (1 - t)) / sum(w),
      treat_all = (tp + fn - (tn + fp) * t / (1 - t)) / sum(w), treat_none = 0)
  }))
}

#' Weighted calibration bins
#' @param y,p,w Outcomes, probabilities and weights.
#' @param bins Number of equal-width probability bins.
#' @return A data frame with count, weight, predicted and observed proportions.
#' @export
svy_calibration_bins <- function(y, p, w = rep(1, length(y)), bins = 10L) {
  y <- .binary(y); p <- .prob(p, length(y)); w <- .positive(w, length(y), TRUE)
  if (length(bins) != 1L || bins < 2L || bins != as.integer(bins)) stop("bins must be an integer >= 2.")
  g <- pmin(bins, floor(p * bins) + 1L)
  do.call(rbind, lapply(split(seq_along(y), g), function(i) {
    sw <- sum(w[i])
    data.frame(bin = g[i[1]], n = length(i), weight = sw,
      predicted = if (sw > 0) sum(w[i] * p[i]) / sw else NA_real_,
      observed = if (sw > 0) sum(w[i] * y[i]) / sw else NA_real_)
  }))
}

#' Weighted empirical risk percentiles
#' @param p Probabilities.
#' @param w Nonnegative reference weights.
#' @return Mid-distribution percentiles between zero and one.
#' @export
svy_percentile <- function(p, w = rep(1, length(p))) {
  p <- .prob(p, length(p)); w <- .positive(w, length(p), TRUE)
  o <- order(p); g <- cumsum(c(TRUE, diff(p[o]) != 0))
  mass <- as.numeric(rowsum(w[o], g, reorder = FALSE))
  ans <- numeric(length(p)); ans[o] <- ((cumsum(mass) - mass / 2) / sum(w))[g]; ans
}
