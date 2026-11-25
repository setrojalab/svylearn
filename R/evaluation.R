#' Evaluate the frozen survey-learning metric suite
#'
#' All 12 core metric keys have stable locations in the returned object. Metrics
#' needing additional data remain NA or NULL with an availability explanation.
#' This does not manufacture transport, sensitivity or explanation evidence from
#' prediction scores alone. A fitted object evaluates its consensus OOF predictions;
#' repeat-average performance remains available in the fitted object's metrics.
#' @param design A svy_learn_design, or a fitted svy_learn_fit.
#' @param y,p Outcomes and predictions. Omit both when passing a fitted object.
#' @param family binomial, gaussian, poisson or multinomial.
#' @param iid_predictions Optional independently generated IID OOF predictions
#'   for the same rows and response schema. The caller establishes their provenance.
#' @param thresholds Optional binary decision thresholds; NULL skips utility.
#' @param population_size Optional known population size for projected counts.
#' @param group Optional domain labels for binary calibration audits.
#' @param shap Optional externally computed local SHAP matrix.
#' @param rake,transport,conformal,sensitivity,importance Optional results from
#'   svy_rake_audit, svy_transport, svy_conformal_eval, svy_sensitivity and
#'   svy_importance_infer, respectively.
#' @param active_predictors Optional prespecified active predictor count for
#'   heuristic support flags; not an automatic model degrees-of-freedom estimate.
#' @param min_events Minimum domain event and nonevent counts.
#' @param level Confidence level.
#' @param replicate_design Optional matching survey replicate design. ICI/E90
#'   replication refits the smoother; no model training is repeated here.
#' @param interpret Whether printing includes the glossary interpretations.
#' @return An object of class svy_eval containing named metrics, inference tables,
#'   availability, support flags, provenance and interpretation metadata.
#' @export
#' @examples
#' d <- svy_example(200)
#' z <- svy_design(d, "weight", "psu", "stratum")
#' if (requireNamespace("survey", quietly = TRUE)) {
#'   audit <- svy_eval(z, d$y, plogis(-2 + d$x1), thresholds = .1)
#'   audit$loss_svy
#' }
svy_eval <- function(design, y = NULL, p = NULL, family = "binomial", iid_predictions = NULL,
                      thresholds = NULL, population_size = NULL, group = NULL, shap = NULL,
                      rake = NULL, transport = NULL, conformal = NULL, sensitivity = NULL,
                      importance = NULL, active_predictors = NULL, min_events = 10L,
                      level = .95, replicate_design = NULL, interpret = FALSE) {
  prediction_source <- "supplied predictions; caller establishes out-of-sample provenance"
  if (inherits(design, "svy_learn_fit")) {
    fit <- design
    if (!is.null(y) || !is.null(p)) stop("Do not override outcomes or predictions when evaluating a fitted object.")
    design <- fit$design; y <- design$data[[fit$outcome]]; p <- fit$consensus
    family <- if (is.null(fit$family)) "binomial" else fit$family
    prediction_source <- "consensus out-of-fold predictions; not the mean of repeat-specific losses"
    if (inherits(fit, "svy_attrition_fit")) prediction_source <- paste(prediction_source,
      "and mean cross-fitted retention weights; fixed-weight inference omits retention-model estimation variance")
  }
  family <- match.arg(family, c("binomial", "gaussian", "poisson", "multinomial")); .eval_level(level)
  if (is.null(y) || is.null(p) || length(y) != nrow(design$data)) stop("Supply outcomes and predictions matching design rows.")
  if (!is.null(active_predictors) && (length(active_predictors) != 1L ||
      !is.finite(active_predictors) || active_predictors < 0 || active_predictors != as.integer(active_predictors)))
    stop("active_predictors must be a nonnegative integer.")
  losses <- .eval_losses(y, p, family); n <- length(y); w <- design$w
  opt <- options(survey.lonely.psu = design$lonely.psu); on.exit(options(opt), add = TRUE)
  des <- .eval_design(design, replicate_design)
  scores <- .eval_mean(des, losses, level); scores$metric[1L] <- "loss_svy"
  raw_variance <- if (n > 1L) stats::var(losses[, 1L]) / n else NA_real_
  deff <- if (is.finite(raw_variance) && raw_variance > 0) scores$se[1L]^2 / raw_variance else NA_real_
  neff <- if (is.finite(deff) && deff > 0) n / deff else if (identical(deff, 0)) Inf else NA_real_
  diagnostics <- svy_diagnostics(design)
  singleton <- list(strata = diagnostics$singleton_strata, rule = design$lonely.psu,
    sampling_df = diagnostics$design_df, variance_df = survey::degf(des),
    variance_engine = if (inherits(des, "svyrep.design")) "supplied replicate conventions" else "survey linearization",
    interpretation = "Audit of specified policy; not a regulatory compliance certification")
  cal <- .eval_calibration(des, y, p, family, level)
  getcal <- function(key) if (is.null(cal)) NA_real_ else cal$estimate[match(key, cal$metric)]
  contrast <- NA_real_
  if (!is.null(iid_predictions)) contrast <- mean(.eval_losses(y, iid_predictions, family)[, 1L]) - scores$estimate[1L]
  point <- svy_metrics(y, p, w, family)
  ici <- e90 <- ipa <- NA_real_; smoother <- NULL; nonlinear <- NULL
  if (family == "binomial") {
    ce <- svy_calibration_error(y, p, w); ici <- ce$ICI; e90 <- ce$E90; smoother <- ce$smoother_status
    ipa <- unname(point["scaled_brier"])
    if (is.finite(ipa)) {
      est <- survey::svymean(cbind(brier = (y - p)^2, prevalence = y), des)
      m <- stats::coef(est); v <- stats::vcov(est); b <- m[1L]; mu <- m[2L]; den <- mu * (1 - mu)
      gradient <- c(-1 / den, b * (1 - 2 * mu) / den^2)
      se <- sqrt(max(0, as.numeric(t(gradient) %*% v %*% gradient)))
      df <- survey::degf(des); q <- if (df > 0) stats::qt((1 + level) / 2, df) else NA_real_
      nonlinear <- data.frame(metric = "ipa_svy", estimate = ipa, se = se,
        ci_lower = ipa - q * se, ci_upper = ipa + q * se, df = df, method = "delta; estimated prevalence reference")
    }
    if (inherits(des, "svyrep.design")) {
      ic <- .eval_replicate(des, function(ww) {
        keep <- ww > 0
        cc <- svy_calibration_error(y[keep], p[keep], ww[keep])
        c(ici_svy = cc$ICI, e90_svy = cc$E90)
      }, level)
      nonlinear <- rbind(nonlinear, ic)
    }
  }
  utility <- if (!is.null(thresholds) && family == "binomial")
    svy_decision_infer(design, y, p, thresholds, population_size, level, replicate_design) else NULL
  domains <- if (!is.null(group) && family == "binomial")
    svy_domain_calibration(design, y, p, group, min_events, level, replicate_design) else NULL
  sh <- if (is.null(shap)) NULL else svy_shap(design, shap, level, replicate_design)
  supplied <- list(rake = rake, transport = transport, conformal = conformal,
    sensitivity = sensitivity, importance = importance)
  fields <- c(rake = "rake_residual", transport = "transport_smd", conformal = "conf_coverage",
    sensitivity = "sensitivity_gamma", importance = "vip_psu")
  for (nm in names(supplied)) if (!is.null(supplied[[nm]]) &&
      (!is.list(supplied[[nm]]) || !fields[[nm]] %in% names(supplied[[nm]])))
    stop("Invalid ", nm, " result; use the documented audit function.")
  val <- function(x, key, default = NA_real_) if (is.null(x) || is.null(x[[key]])) default else x[[key]]
  decision <- function(key) if (is.null(utility)) NULL else utility[utility$metric == key, , drop = FALSE]
  flags <- character()
  if (!is.null(active_predictors) && is.finite(neff) && neff < 10 * active_predictors)
    flags <- c(flags, "Heuristic: predictive effective n < 10 * active_predictors")
  if (family == "binomial" && is.finite(neff) && neff * sum(w * y) / sum(w) < 20)
    flags <- c(flags, "Heuristic: predictive effective n times prevalence < 20")
  availability <- svy_glossary()[c("core", "name", "key", "required_inputs")]
  available <- c(TRUE, is.finite(ipa), !is.null(iid_predictions), is.finite(deff), !is.na(neff), TRUE,
    !is.null(cal), is.finite(ici), !is.null(rake), !is.null(transport), TRUE, !is.null(conformal),
    !is.null(sensitivity), !is.null(utility), !is.null(utility), !is.null(domains), !is.null(sh) || !is.null(importance))
  availability$status <- ifelse(available, "computed", "not computed: see required_inputs")
  if (xor(!is.null(sh), !is.null(importance))) availability$status[17L] <- "partially computed: SHAP and VIP require separate inputs"
  if (!is.null(domains) && is.null(domains$coefficients)) availability$status[16L] <- "not estimable: inadequate domain support"
  if (family != "binomial") availability$status[c(2, 8, 14, 15, 16)] <- "not applicable to this outcome family"
  result <- list(loss_svy = scores$estimate[1L], ipa_svy = ipa, optimism_bias = contrast,
    deff_loss = deff, n_eff_pred = neff, singleton_audit = singleton,
    cal_intercept = getcal("cal_intercept"), cal_joint_intercept = getcal("cal_joint_intercept"),
    cal_slope = getcal("cal_slope"), ici_svy = ici, e90_svy = e90,
    rake_residual = val(rake, "rake_residual"), transport_smd = val(transport, "transport_smd", NULL),
    transport_div = val(transport, "transport_div"), conf_coverage = val(conformal, "conf_coverage"),
    conf_width = val(conformal, "conf_width"), sensitivity_gamma = val(sensitivity, "sensitivity_gamma"),
    net_benefit_svy = decision("net_benefit_svy"), pop_tp = decision("pop_tp"), pop_fp = decision("pop_fp"),
    pop_tn = decision("pop_tn"), pop_fn = decision("pop_fn"), nns_svy = decision("nns_svy"),
    nns_population = decision("nns_population"), domain_cal_gap = val(domains, "domain_cal_gap", NULL),
    shap_svy = sh, vip_psu = val(importance, "vip_psu", NULL),
    ci_lower = scores$ci_lower[1L], ci_upper = scores$ci_upper[1L],
    intervals = rbind(scores, cal, nonlinear), utility = utility, calibration = cal,
    domains = domains, diagnostics = diagnostics, availability = availability, support_flags = flags,
    point_metrics = point, supplied_audits = supplied, smoother_status = smoother,
    family = family, prediction_source = prediction_source, interpret = isTRUE(interpret),
    variance_reference = "Unweighted sample variance of row loss divided by nominal n; SRS with replacement",
    uncertainty_target = "Conditional on supplied predictions; no training/tuning refits; pointwise intervals unless stated")
  structure(result, class = "svy_eval")
}

#' @rdname svy_eval
#' @param ... Arguments passed to svy_eval or to an S3 method.
#' @export
evaluate <- function(...) svy_eval(...)

#' @export
print.svy_eval <- function(x, interpret = x$interpret, ...) {
  cat("Survey evaluation (", x$family, ")\n", sep = "")
  print(x$intervals, row.names = FALSE)
  cat("Loss DEFF:", x$deff_loss, " Predictive effective n:", x$n_eff_pred, "\n")
  cat(x$prediction_source, "\n", x$uncertainty_target, "\n")
  if (length(x$support_flags)) cat(paste(x$support_flags, collapse = "\n"), "\n")
  if (isTRUE(interpret)) {
    g <- svy_glossary()
    for (i in seq_len(nrow(g))) cat(g$name[i], ": ", g$interpretation[i],
      " [", x$availability$status[i], "]\n", sep = "")
  }
  invisible(x)
}

#' @export
summary.svy_eval <- function(object, interpret = object$interpret, ...) {
  object$interpret <- isTRUE(interpret); object
}
