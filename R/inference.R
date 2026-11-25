#' Design-based evaluation of fixed predictions
#'
#' Uncertainty is conditional on the supplied prediction vector. It does not
#' include training, tuning or recalibration uncertainty. Use a full repeated
#' sampling experiment to assess unconditional pipeline coverage. Domains use
#' the parent design, retaining its sampling structure.
#' @param design A svy_learn_design.
#' @param y,p Outcomes and predictions in data row order, as in svy_metrics.
#' @param domain Optional logical domain indicator.
#' @param level Confidence level.
#' @param replicate_design Optional survey replicate design in the same row order.
#' @param family binomial, gaussian, poisson or multinomial. Nonbinary families
#'   return design-based mean-loss intervals; supply replicate weights within
#'   design for replicate-based mean-loss inference.
#' @return A list with proper-score inference, calibration regression, point
#'   metrics, optional rank-metric replicate inference, and the uncertainty target.
#' @export
svy_infer <- function(design, y, p, domain = rep(TRUE, length(y)), level = .95,
                      replicate_design = NULL, family = "binomial") {
  family <- match.arg(family, c("binomial", "gaussian", "poisson", "multinomial"))
  if (family != "binomial") return(.response_infer(design, y, p, domain, level,
    replicate_design, family))
  .need("survey"); y <- .binary(y); p <- .prob(p, length(y))
  if (length(y) != nrow(design$data)) stop("Row count mismatch.")
  if (!is.logical(domain) || length(domain) != length(y) || anyNA(domain) || !any(domain)) stop("Invalid domain.")
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) stop("Invalid confidence level.")
  opt <- options(survey.lonely.psu = design$lonely.psu)
  on.exit(options(opt), add = TRUE)
  des <- as_survey_design(design); pc <- .clip(p)
  des$variables$.svyl_logloss <- -(y * log(pc) + (1 - y) * log1p(-pc))
  des$variables$.svyl_brier <- (y - p)^2
  des$variables$.svyl_y <- y; des$variables$.svyl_lp <- stats::qlogis(pc)
  des <- .survey_domain(des, domain)
  df <- survey::degf(des)
  if (df <= 0) stop("Positive survey degrees of freedom are required for intervals.")
  est <- survey::svymean(~.svyl_logloss + .svyl_brier, des, deff = "replace")
  val <- stats::coef(est); se <- sqrt(diag(stats::vcov(est)))
  q <- stats::qt((1 + level) / 2, df)
  scores <- data.frame(metric = c("log_loss", "brier"), estimate = as.numeric(val),
    se = se, lower = val - q * se, upper = val + q * se,
    deff = as.numeric(survey::deff(est)), df = df, row.names = NULL)
  cal <- NULL
  if (length(unique(y[domain])) == 2L && length(unique(p[domain])) > 1L) {
    fit <- survey::svyglm(.svyl_y ~ .svyl_lp, design = des, family = stats::quasibinomial())
    citl <- survey::svyglm(.svyl_y ~ 1 + offset(.svyl_lp),
                          design = des, family = stats::quasibinomial())
    cval <- c(stats::coef(citl)[1], stats::coef(fit)[2])
    cse <- c(sqrt(stats::vcov(citl)[1, 1]), sqrt(stats::vcov(fit)[2, 2]))
    cdf <- c(citl$df.residual, fit$df.residual)
    cq <- stats::qt((1 + level) / 2, pmax(cdf, 1))
    cal <- data.frame(parameter = c("calibration_in_the_large", "calibration_slope"),
      estimate = as.numeric(cval), se = cse, lower = cval - cq * cse,
      upper = cval + cq * cse, df = cdf, row.names = NULL)
  }
  ranks <- NULL
  if (!is.null(replicate_design)) {
    rd <- replicate_design
    if (nrow(rd$variables) != length(y)) stop("Replicate-design row count mismatch.")
    rd$variables$.svyl_y <- y; rd$variables$.svyl_p <- p; rd <- .survey_domain(rd, domain)
    rr <- survey::withReplicates(rd, function(w, data)
      svy_metrics(data$.svyl_y, data$.svyl_p, w)[c("auc", "average_precision")])
    rv <- stats::coef(rr); rs <- sqrt(diag(stats::vcov(rr)))
    rq <- stats::qt((1 + level) / 2, survey::degf(rd))
    ranks <- data.frame(metric = names(rv), estimate = as.numeric(rv), se = rs,
                        lower = rv - rq * rs, upper = rv + rq * rs, row.names = NULL)
  }
  list(scores = scores, calibration = cal, ranks = ranks,
    metrics = svy_metrics(y[domain], p[domain], design$w[domain]),
    target = "Design-based uncertainty conditional on supplied predictions; no model refitting")
}

#' Domain performance audit with support suppression
#' @param design A svy_learn_design.
#' @param y,p Outcomes and probabilities.
#' @param group Group labels in data row order.
#' @param min_events Minimum events and nonevents for binary calibration inference;
#'   minimum count per declared class for multinomial inference, or minimum
#'   domain sample size for gaussian and poisson inference.
#' @param family binomial, gaussian, poisson or multinomial.
#' @return A list of domain-specific inference and support information.
#' @export
svy_subgroups <- function(design, y, p, group, min_events = 10L, family = "binomial") {
  family <- match.arg(family, c("binomial", "gaussian", "poisson", "multinomial"))
  if (family != "binomial") {
    y <- .response_y(y, family); p <- .response_p(p, length(y), family, levels(y))
    if (length(y) != nrow(design$data) || length(group) != length(y) || anyNA(group)) stop("Invalid groups or row count.")
    if (length(min_events) != 1L || !is.finite(min_events) || min_events < 1 || min_events != as.integer(min_events)) stop("Invalid support threshold.")
    return(lapply(split(seq_along(y), group), function(i) {
      counts <- if (is.factor(y)) table(y[i]) else NULL
      sufficient <- if (is.null(counts)) length(i) >= min_events else all(counts >= min_events)
      domain <- seq_along(y) %in% i
      list(n = length(i), class_counts = counts, sufficient = sufficient,
        metrics = svy_metrics(y[i], if (is.matrix(p)) p[i, , drop = FALSE] else p[i], design$w[i], family),
        inference = if (sufficient) svy_infer(design, y, p, domain = domain, family = family) else NULL)
    }))
  }
  y <- .binary(y); p <- .prob(p, length(y))
  if (length(group) != length(y) || anyNA(group)) stop("Invalid groups.")
  lapply(split(seq_along(y), group), function(i) {
    sufficient <- sum(y[i]) >= min_events && sum(1 - y[i]) >= min_events
    domain <- seq_along(y) %in% i
    list(n = length(i), events = sum(y[i]), sufficient = sufficient,
      metrics = svy_metrics(y[i], p[i], design$w[i]),
      inference = if (sufficient) svy_infer(design, y, p, domain = domain) else NULL)
  })
}

#' Aggregate provided SHAP contributions with survey weights
#' @param design A svy_learn_design.
#' @param shap Numeric matrix of externally computed SHAP contributions.
#' @param level Confidence level.
#' @param replicate_design Optional matching survey replicate design.
#' @return Design-based mean absolute attributions, standard errors and intervals.
#'   Their interpretation depends on the supplied model and SHAP background.
#' @export
svy_shap <- function(design, shap, level = .95, replicate_design = NULL) {
  .need("survey")
  if (!is.matrix(shap) || !is.numeric(shap) || nrow(shap) != nrow(design$data) ||
      any(!is.finite(shap))) stop("Invalid SHAP matrix.")
  opt <- options(survey.lonely.psu = design$lonely.psu); on.exit(options(opt), add = TRUE)
  .eval_level(level)
  des <- .eval_design(design, replicate_design)
  est <- .eval_mean(des, abs(shap), level)
  data.frame(feature = if (is.null(colnames(shap))) paste0("x", seq_len(ncol(shap))) else colnames(shap),
    mean_absolute_shap = est$estimate, shap_svy = est$estimate,
    se = est$se, ci_lower = est$ci_lower, ci_upper = est$ci_upper,
    df = est$df, method = est$method, row.names = NULL)
}

#' Out-of-fold permutation importance within sampling clusters
#' @param fit A svy_learn_fit.
#' @param repeats Number of permutation draws.
#' @param seed Random seed.
#' @return Data frame of primary-loss increases (log loss for classification,
#'   squared error for gaussian, Poisson deviance for counts). These are associational
#'   perturbation summaries, not causal effects or conditional-independence tests.
#' @export
svy_importance <- function(fit, repeats = 5L, seed = 1L) {
  .seed(seed, {
    d <- fit$design; y <- d$data[[fit$outcome]]; out <- list()
    generalized <- inherits(fit, "svy_response_fit")
    family <- if (generalized) fit$family else "binomial"
    metric <- colnames(fit$metrics)[1L]
    for (r in seq_len(ncol(fit$folds))) for (j in seq_len(repeats)) for (variable in fit$predictors) {
      nd <- d$data
      for (i in split(seq_len(nrow(nd)), d$cluster)) nd[[variable]][i] <- nd[[variable]][i[sample.int(length(i))]]
      p <- if (family == "multinomial") matrix(0, nrow(nd), length(fit$levels),
        dimnames = list(NULL, fit$levels)) else numeric(nrow(nd))
      for (k in seq_len(max(fit$folds))) {
        i <- which(fit$folds[, r] == k)
        f <- fit$outer_fits[[(r - 1L) * max(fit$folds) + k]]
        pp <- if (generalized) .response_pipeline_predict(f, nd[i, , drop = FALSE]) else
          .pipeline_predict(f, nd[i, , drop = FALSE])
        if (family == "multinomial") p[i, ] <- pp else p[i] <- pp
      }
      delta <- svy_metrics(y, p, d$w, family)[metric] - fit$metrics[r, metric]
      out[[length(out) + 1L]] <- data.frame(variable = variable, repeat_id = r,
        permutation = j, log_loss_increase = as.numeric(delta))
      if (family %in% c("gaussian", "poisson")) names(out[[length(out)]])[4L] <- paste0(metric, "_increase")
    }
    do.call(rbind, out)
  })
}
