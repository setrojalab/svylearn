#' Joint domain calibration and maximal departures from ideal calibration
#' @param design A svy_learn_design.
#' @param y,p Binary outcomes and probabilities.
#' @param group Nonmissing domain labels in row order.
#' @param min_events Minimum events and nonevents in each included domain.
#' @param level Confidence level.
#' @param replicate_design Optional matching survey replicate design.
#' @return Domain coefficients, joint covariance, pairwise contrasts and
#'   domain_cal_gap (maximum absolute joint intercept and slope departure from one).
#'   Maximal-gap intervals use a Bonferroni coefficient rectangle; pairwise
#'   intervals use the joint covariance and are pointwise. These are descriptive
#'   calibration audits, not proof of algorithmic fairness or causal parity.
#' @export
svy_domain_calibration <- function(design, y, p, group, min_events = 10L,
                                    level = .95, replicate_design = NULL) {
  y <- .binary(y); p <- .prob(p, length(y)); .eval_level(level)
  if (length(y) != nrow(design$data) || length(group) != length(y) || anyNA(group))
    stop("Domain labels and outcomes must match the design rows.")
  if (length(min_events) != 1L || !is.finite(min_events) || min_events < 1 ||
      min_events != as.integer(min_events)) stop("Invalid minimum event count.")
  group <- as.character(group); labels <- sort(unique(group))
  support <- do.call(rbind, lapply(labels, function(g) {
    i <- group == g
    data.frame(domain = g, n = sum(i), events = sum(y[i]), nonevents = sum(1 - y[i]),
      sufficient = min(sum(y[i]), sum(1 - y[i])) >= min_events && length(unique(p[i])) > 1L)
  }))
  labels <- support$domain[support$sufficient]
  if (!length(labels)) return(list(domain_cal_gap = c(intercept = NA_real_, slope = NA_real_),
    support = support, coefficients = NULL, contrasts = NULL, covariance = NULL,
    intervals = NULL, status = "No domain has sufficient support"))
  opt <- options(survey.lonely.psu = design$lonely.psu); on.exit(options(opt), add = TRUE)
  des <- .eval_design(design, replicate_design)
  des$variables$.svye_y <- y; des$variables$.svye_lp <- stats::qlogis(.clip(p))
  des$variables$.svye_group <- factor(group, levels = labels)
  des <- .survey_domain(des, group %in% labels)
  if (length(labels) == 1L) {
    fit <- survey::svyglm(.svye_y ~ .svye_lp, des, family = stats::quasibinomial())
    L <- diag(2)
  } else {
    fit <- survey::svyglm(.svye_y ~ 0 + .svye_group + .svye_group:.svye_lp,
      des, family = stats::quasibinomial())
    nd <- data.frame(.svye_group = factor(labels, levels = labels), .svye_lp = 0)
    form <- stats::delete.response(stats::terms(fit))
    z0 <- stats::model.matrix(form, nd, contrasts.arg = fit$contrasts)
    nd$.svye_lp <- 1; z1 <- stats::model.matrix(form, nd, contrasts.arg = fit$contrasts)
    L <- matrix(0, 2L * length(labels), ncol(z0))
    L[seq(1, nrow(L), 2), ] <- z0; L[seq(2, nrow(L), 2), ] <- z1 - z0
  }
  if (!isTRUE(fit$converged) || any(!is.finite(stats::coef(fit))))
    stop("Joint domain calibration is not estimable; increase support or revise domains.")
  beta <- as.numeric(L %*% stats::coef(fit)); V <- L %*% stats::vcov(fit) %*% t(L)
  se <- sqrt(pmax(0, diag(V))); df <- fit$df.residual
  q <- if (df > 0) stats::qt((1 + level) / 2, df) else NA_real_
  coefficients <- data.frame(domain = rep(labels, each = 2L),
    parameter = rep(c("joint_intercept", "slope"), length(labels)),
    estimate = beta, se = se, ci_lower = beta - q * se, ci_upper = beta + q * se, df = df)
  # One joint fit retains cross-domain covariance, including domains sharing PSUs.
  contrasts <- list()
  if (length(labels) > 1L) for (a in seq_len(length(labels) - 1L))
    for (b in seq.int(a + 1L, length(labels))) for (j in 1:2) {
      ia <- 2L * (a - 1L) + j; ib <- 2L * (b - 1L) + j
      ss <- sqrt(max(0, V[ia, ia] + V[ib, ib] - 2 * V[ia, ib])); value <- beta[ia] - beta[ib]
      contrasts[[length(contrasts) + 1L]] <- data.frame(domain_a = labels[a], domain_b = labels[b],
        parameter = c("joint_intercept", "slope")[j], estimate = value, se = ss,
        ci_lower = value - q * ss, ci_upper = value + q * ss, df = df)
    }
  difference <- beta - rep(c(0, 1), length(labels))
  qb <- if (df > 0) stats::qt(1 - (1 - level) / (2 * length(beta)), df) else NA_real_
  low <- difference - qb * se; high <- difference + qb * se
  gaps <- vapply(1:2, function(j) max(abs(difference[seq(j, length(beta), 2)])), numeric(1))
  names(gaps) <- c("intercept", "slope")
  intervals <- do.call(rbind, lapply(1:2, function(j) {
    i <- seq(j, length(beta), 2)
    data.frame(parameter = names(gaps)[j], estimate = gaps[j],
      ci_lower = max(ifelse(low[i] <= 0 & high[i] >= 0, 0, pmin(abs(low[i]), abs(high[i])))),
      ci_upper = max(abs(c(low[i], high[i]))), method = "Bonferroni coefficient rectangle")
  }))
  list(domain_cal_gap = gaps, support = support, coefficients = coefficients,
    covariance = V, contrasts = if (length(contrasts)) do.call(rbind, contrasts) else NULL,
    intervals = intervals, status = "Joint-model calibration; conditional predictions")
}
