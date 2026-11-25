#' Weighted ROC and precision-recall curves
#' @param y,p,w Binary outcomes, probabilities and weights.
#' @param max_fpr Maximum false-positive rate for the unstandardized partial AUC.
#' @return A list with threshold coordinates and partial AUC. Both outcome
#'   classes must have positive total weight. Tied scores enter together.
#' @export
svy_curves <- function(y, p, w = rep(1, length(y)), max_fpr = .1) {
  y <- .binary(y); p <- .prob(p, length(y)); w <- .positive(w, length(y), TRUE)
  if (length(max_fpr) != 1L || !is.finite(max_fpr) || max_fpr <= 0 || max_fpr > 1)
    stop("max_fpr must be in (0, 1].")
  posw <- sum(w * y); negw <- sum(w * (1 - y))
  if (posw <= 0 || negw <= 0) stop("Both outcome classes need positive total weight.")
  o <- order(p, decreasing = TRUE); g <- cumsum(c(TRUE, diff(p[o]) != 0))
  pos <- as.numeric(rowsum(w[o] * y[o], g, reorder = FALSE))
  neg <- as.numeric(rowsum(w[o] * (1 - y[o]), g, reorder = FALSE))
  tp <- c(0, cumsum(pos)); fp <- c(0, cumsum(neg))
  precision <- ifelse(tp + fp > 0, tp / (tp + fp), 1)
  curve <- data.frame(threshold = c(Inf, unique(p[o])), fpr = fp / negw,
    tpr = tp / posw, recall = tp / posw, precision = precision)
  # Integrate linear ROC segments, truncating the crossing segment exactly.
  area <- 0
  for (i in seq_len(nrow(curve) - 1L)) {
    a <- curve$fpr[i]; b <- curve$fpr[i + 1L]
    if (a >= max_fpr || b <= a) next
    right <- min(max_fpr, b)
    t <- curve$tpr[i] + (curve$tpr[i + 1L] - curve$tpr[i]) * (right - a) / (b - a)
    area <- area + (right - a) * (curve$tpr[i] + t) / 2
  }
  list(curve = curve, partial_auc = area, max_fpr = max_fpr)
}

#' Weighted smooth calibration error diagnostics
#' @param y,p,w Outcomes, probabilities and weights from an evaluation sample.
#' @param df Degrees of freedom of the natural spline calibration smoother.
#' @return Named ICI, E50 and E90 diagnostics, and the smoother status. These
#'   depend on the smoother; they are not confidence bounds or corrected bias.
#' @export
svy_calibration_error <- function(y, p, w = rep(1, length(y)), df = 3L) {
  y <- .binary(y); p <- .prob(p, length(y)); w <- .positive(w, length(y))
  f <- svy_calibrate(y, p, w, method = "spline", df = df)
  error <- abs(stats::predict(f, p) - p)
  list(ICI = sum(w * error) / sum(w), E50 = .wquantile(error, w, .5),
       E90 = .wquantile(error, w, .9), smoother_status = f$status)
}
