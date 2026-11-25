#' Sharp bounded-weight mean sensitivity analysis
#'
#' The uncertainty set permits each analysis weight to be multiplied by any
#' number in [1/gamma, gamma]. Sorting values reduces minimization of the weighted
#' mean to a threshold tilt. The reported point bounds are sharp for this box.
#' This is not an unmeasured-confounding E-value or a missing-outcome model.
#' @param values Finite row-level utility, loss, or calibration residual values.
#' @param w Nonnegative baseline weights.
#' @param gamma Sensitivity multipliers at least one.
#' @return Lower and upper bounds and gamma. No sampling intervals are implied.
#' @export
svy_weight_bounds <- function(values, w = rep(1, length(values)), gamma = c(1, 1.5, 2, 3, 5)) {
  if (!is.numeric(values) || !length(values) || any(!is.finite(values))) stop("values must be finite numeric values.")
  w <- .positive(w, length(values), TRUE)
  if (!length(gamma) || any(!is.finite(gamma)) || any(gamma < 1)) stop("gamma must be finite and at least one.")
  calc <- .tilt_setup(values, w)
  data.frame(gamma = gamma, lower = vapply(gamma, function(g) calc(g)$estimate, numeric(1)),
    upper = vapply(gamma, function(g) - .tilt_setup(-values, w)(g)$estimate, numeric(1)))
}
.tilt_setup <- function(values, w) {
  o <- order(values); v <- values[o]; ww <- w[o] / sum(w)
  cw <- c(0, cumsum(ww)); cv <- c(0, cumsum(ww * v)); total <- sum(ww * v)
  function(gamma, multiplier = FALSE) {
    # At an extremum all values below the resulting mean get the upper multiplier.
    # Enumerating n+1 thresholds avoids numerical optimization and proves sharpness.
    high <- gamma; low <- 1 / gamma
    means <- (low * total + (high - low) * cv) / (low + (high - low) * cw)
    k <- which.min(means) - 1L
    ans <- list(estimate = means[k + 1L])
    if (multiplier) {
      u <- rep(low, length(v)); if (k > 0) u[seq_len(k)] <- high
      ans$multiplier <- u[order(o)]
    }
    ans
  }
}

#' Tipping analysis for bounded survey-weight distortion
#' @param design A svy_learn_design.
#' @param values Finite row-level utility (for example binary net-benefit contributions).
#' @param gamma Increasing search grid starting at one.
#' @param null_value Value at which utility is considered overturned.
#' @param target point_bound or tilted_lower_ci. The latter follows the exact
#'   mean-minimizing tilt and recomputes a conditional interval there; it does
#'   not globally minimize the lower confidence bound over all possible tilts.
#' @param level Confidence level for tilted_lower_ci.
#' @param tolerance Numerical root tolerance.
#' @return sensitivity_gamma, search curve, bracket and target. Infinity means
#'   no crossing was detected in the supplied range. For the CI path, the first
#'   grid-bracketed crossing is reported; global monotonicity is not assumed.
#' @export
svy_sensitivity <- function(design, values, gamma = exp(seq(0, log(5), length.out = 41)),
                            null_value = 0, target = "point_bound", level = .95, tolerance = 1e-8) {
  target <- match.arg(target, c("point_bound", "tilted_lower_ci")); .eval_level(level)
  if (!is.numeric(values) || length(values) != nrow(design$data) || any(!is.finite(values)))
    stop("values must be finite and match the design rows.")
  if (length(gamma) < 2L || any(!is.finite(gamma)) || gamma[1L] != 1 || any(diff(gamma) <= 0))
    stop("gamma must be an increasing grid starting at one.")
  if (length(null_value) != 1L || !is.finite(null_value) || length(tolerance) != 1L ||
      !is.finite(tolerance) || tolerance <= 0) stop("Invalid null value or tolerance.")
  calc <- .tilt_setup(values, design$w)
  opt <- options(survey.lonely.psu = design$lonely.psu); on.exit(options(opt), add = TRUE)
  if (target == "tilted_lower_ci") des <- .eval_design(design)
  at <- function(g) {
    tilt <- calc(g, multiplier = target == "tilted_lower_ci")
    if (target == "point_bound") return(c(bound = tilt$estimate, estimate = tilt$estimate, se = NA_real_))
    # Ratio of tilted totals retains the original sampling design and its FPCs.
    # The tilt is held fixed for this conditional interval.
    stat <- survey::svyratio(matrix(values * tilt$multiplier, ncol = 1L),
      matrix(tilt$multiplier, ncol = 1L), des)
    tab <- .eval_table(stat, survey::degf(des), level, "fixed adversarial tilt")
    if (!is.finite(tab$ci_lower)) stop("The tilted confidence bound requires positive design degrees of freedom.")
    c(bound = tab$ci_lower, estimate = tab$estimate, se = tab$se)
  }
  curve <- data.frame(gamma = gamma, t(vapply(gamma, at, numeric(3))))
  hit <- which(curve$bound <= null_value)
  bracket <- c(NA_real_, NA_real_); tipping <- Inf
  if (length(hit)) {
    j <- hit[1L]
    if (j == 1L) { tipping <- 1; bracket <- c(1, 1) } else {
      bracket <- gamma[c(j - 1L, j)]
      tipping <- stats::uniroot(function(g) at(g)[["bound"]] - null_value,
        bracket, tol = tolerance)$root
    }
  }
  list(sensitivity_gamma = tipping, target = target, curve = curve, bracket = bracket,
    searched_range = range(gamma), global_point_optimum = target == "point_bound",
    interpretation = if (target == "point_bound") "Sharp box-weight mean bound; no sampling CI" else
      "First grid-bracketed crossing along the mean-minimizing tilt; conditional CI, not globally optimal CI tipping")
}
