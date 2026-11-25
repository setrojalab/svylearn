#' Fit a weighted probability recalibration map
#' @param y,p,w Calibration outcomes, honest predictions and weights.
#' @param method platt, isotonic, spline or none.
#' @param df Natural spline degrees of freedom.
#' @param min_events Minimum raw and Kish effective counts in each class for
#'   estimating a calibration slope. Below this, use a smoothed intercept shift.
#' @return A fitted svy_calibrator with a predict method. Single-class training
#'   samples use the weighted event proportion and record a fallback status.
#' @export
svy_calibrate <- function(y, p, w = rep(1, length(y)), method = "platt", df = 3L,
                          min_events = 5L) {
  y <- .binary(y); p <- .prob(p, length(y)); w <- .positive(w, length(y))
  method <- match.arg(method, c("platt", "isotonic", "spline", "none"))
  out <- list(method = method, status = "ok", df = df)
  if (method == "none") return(structure(out, class = "svy_calibrator"))
  neff <- sum(w)^2 / sum(w^2); prev <- sum(w * y) / sum(w)
  smoothed <- (neff * prev + .5) / (neff + 1)
  effective <- vapply(0:1, function(cls) {
    ww <- w[y == cls]; if (length(ww)) sum(ww)^2 / sum(ww^2) else 0
  }, numeric(1))
  if (min(c(sum(y), sum(1 - y), effective)) < min_events) {
    shifted <- svy_align_prevalence(p, w, smoothed)
    out$offset <- attr(shifted, "offset")
    out$status <- "sparse calibration smoothed intercept-only fallback"
    return(structure(out, class = "svy_calibrator"))
  }
  if (length(unique(y)) < 2L || length(unique(p)) < 2L) {
    out$constant <- smoothed; out$status <- "constant calibration fallback"
  } else if (method == "isotonic") {
    o <- order(p); pp <- p[o]; g <- cumsum(c(TRUE, diff(pp) != 0))
    ww <- as.numeric(rowsum(w[o], g, reorder = FALSE))
    yy <- as.numeric(rowsum(w[o] * y[o], g, reorder = FALSE))
    xx <- unique(pp); starts <- ends <- numeric(length(xx)); val <- mass <- numeric(length(xx)); top <- 0L
    for (j in seq_along(xx)) {
      top <- top + 1L; starts[top] <- ends[top] <- j; val[top] <- yy[j] / ww[j]; mass[top] <- ww[j]
      while (top > 1L && val[top - 1L] > val[top]) {
        v <- (val[top - 1L] * mass[top - 1L] + val[top] * mass[top]) / (mass[top - 1L] + mass[top])
        mass[top - 1L] <- mass[top - 1L] + mass[top]; val[top - 1L] <- v
        ends[top - 1L] <- ends[top]; top <- top - 1L
      }
    }
    fit <- numeric(length(xx))
    for (j in seq_len(top)) fit[starts[j]:ends[j]] <- val[j]
    out$x <- xx; out$fitted <- fit
  } else {
    z <- stats::qlogis(.clip(p)); dd <- data.frame(y = y, z = z, w = w / mean(w))
    if (method == "spline") {
      b <- splines::ns(z, df = df); out$basis <- b; xx <- cbind(1, b)
    } else xx <- cbind(1, z)
    fit <- stats::glm.fit(xx, y, weights = dd$w, family = stats::quasibinomial())
    if (!isTRUE(fit$converged) || any(!is.finite(fit$coefficients))) {
      out$constant <- smoothed; out$status <- "failed calibration fallback"
    } else out$coefficients <- fit$coefficients
  }
  structure(out, class = "svy_calibrator")
}

#' @export
predict.svy_calibrator <- function(object, newdata, ...) {
  p <- .prob(newdata, length(newdata))
  if (object$method == "none") return(p)
  if (!is.null(object$offset)) return(.clip(stats::plogis(stats::qlogis(.clip(p)) + object$offset)))
  if (!is.null(object$constant)) return(rep(.clip(object$constant), length(p)))
  if (object$method == "isotonic") {
    i <- pmax(1L, findInterval(p, object$x)); return(.clip(object$fitted[i]))
  }
  z <- stats::qlogis(.clip(p))
  x <- if (object$method == "spline") cbind(1, stats::predict(object$basis, z)) else cbind(1, z)
  .clip(stats::plogis(as.numeric(x %*% object$coefficients)))
}

#' Align the weighted mean probability with an external prevalence
#' @param p,w Predicted probabilities and target-population weights.
#' @param prevalence Externally known population prevalence in (0, 1).
#' @return Adjusted probabilities, with the fitted log-odds offset as an attribute.
#' @export
svy_align_prevalence <- function(p, w, prevalence) {
  p <- .prob(p, length(p)); w <- .positive(w, length(p))
  if (length(prevalence) != 1L || !is.finite(prevalence) || prevalence <= 0 || prevalence >= 1)
    stop("prevalence must be strictly between zero and one.")
  z <- stats::qlogis(.clip(p))
  a <- stats::uniroot(function(a) sum(w * stats::plogis(z + a)) / sum(w) - prevalence,
                      c(-80, 80), tol = 1e-12)$root
  structure(stats::plogis(z + a), offset = a)
}

#' Calibrate survey weights to known population totals
#' @param design A svy_learn_design.
#' @param formula,population Calibration formula and totals for survey::calibrate.
#' @param ... Additional survey::calibrate arguments.
#' @return A calibrated survey design. This changes weights, not probabilities.
#' @export
svy_rake <- function(design, formula, population, ...) {
  .need("survey")
  survey::calibrate(as_survey_design(design), formula = formula, population = population, ...)
}
