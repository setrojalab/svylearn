#' Audit predicted event totals and calibrated covariate margins
#' @param design A svy_learn_design.
#' @param p Probabilities in design row order.
#' @param known_event_total Externally known positive event total.
#' @param raked_design Optional calibrated survey design returned by svy_rake.
#' @param formula,population Optional margin formula and named target totals.
#' @param level Confidence level.
#' @return rake_residual, event-total inference and optional margin residuals.
#'   Demographic raking alone does not force the event residual to zero.
#' @export
svy_rake_audit <- function(design, p, known_event_total, raked_design = NULL,
                           formula = NULL, population = NULL, level = .95) {
  p <- .prob(p, nrow(design$data)); .eval_level(level)
  if (length(known_event_total) != 1L || !is.finite(known_event_total) || known_event_total <= 0)
    stop("known_event_total must be a positive external total.")
  opt <- options(survey.lonely.psu = design$lonely.psu); on.exit(options(opt), add = TRUE)
  des <- if (is.null(raked_design)) .eval_design(design) else raked_design
  if (!inherits(des, c("survey.design", "svyrep.design")) || nrow(des$variables) != length(p))
    stop("raked_design must be a matching survey design.")
  keys <- unique(c(design$ids, design$strata, design$subject, design$wave))
  if (length(keys) && (!all(keys %in% names(des$variables)) ||
      !all(vapply(keys, function(k) identical(as.character(des$variables[[k]]),
        as.character(design$data[[k]])), logical(1))))) stop("Raked design identifiers do not match.")
  est <- .eval_mean(des, cbind(predicted_event_total = p), level, total = TRUE)
  residual <- est; residual$metric <- "rake_residual"
  for (k in c("estimate", "ci_lower", "ci_upper")) residual[[k]] <- residual[[k]] / known_event_total - 1
  residual$se <- residual$se / known_event_total
  margins <- NULL
  if (!is.null(formula) || !is.null(population)) {
    if (is.null(formula) || is.null(population) || is.null(names(population)) ||
        anyDuplicated(names(population)) || any(!is.finite(population))) stop("Supply a formula and named finite margin totals.")
    x <- stats::model.matrix(formula, des$variables)
    if (!setequal(colnames(x), names(population))) stop("Margin names do not match the model matrix.")
    ww <- if (inherits(des, "svyrep.design")) stats::weights(des, type = "sampling") else stats::weights(des)
    actual <- colSums(x * as.numeric(ww)); target <- population[colnames(x)]
    margins <- data.frame(margin = colnames(x), target = as.numeric(target), actual = actual,
      absolute_residual = actual - target,
      relative_residual = ifelse(target != 0, (actual - target) / target, NA_real_), row.names = NULL)
  }
  list(rake_residual = residual$estimate, event_total = est, intervals = residual, margins = margins,
    target = "Residual against a fixed external event total; supplied calibrated weights and predictions")
}

#' Compare source and target covariate balance and prediction divergence
#' @param source,target Finite numeric matrices with identically named covariates.
#'   Encode categories using a prespecified common indicator schema.
#' @param source_weights,target_weights Nonnegative source and target weights.
#' @param transport_weights Optional final source weights after transport weighting.
#' @param source_prediction,target_prediction Optional prediction vectors, both
#'   evaluated on the same target rows. The first is the transported source rule.
#' @param target_outcome Optional binary target outcomes, requiring target_prediction.
#' @return transport_smd (before/after), transport_div, separate target Brier score
#'   and weight diagnostics. The SMD denominator is the square root of the average
#'   baseline weighted population variances, kept fixed after transport. Neither
#'   covariate balance nor prediction agreement establishes causal transportability.
#' @export
svy_transport <- function(source, target, source_weights = rep(1, nrow(source)),
                           target_weights = rep(1, nrow(target)), transport_weights = NULL,
                           source_prediction = NULL, target_prediction = NULL, target_outcome = NULL) {
  valid <- function(x) is.matrix(x) && is.numeric(x) && nrow(x) > 0 && ncol(x) > 0 &&
    !is.null(colnames(x)) && !anyDuplicated(colnames(x)) && all(is.finite(x))
  if (!valid(source) || !valid(target) || !setequal(colnames(source), colnames(target)))
    stop("Supply finite numeric matrices with matching named covariates.")
  target <- target[, colnames(source), drop = FALSE]
  sw <- .positive(source_weights, nrow(source), TRUE); tw <- .positive(target_weights, nrow(target), TRUE)
  aw <- if (is.null(transport_weights)) sw else .positive(transport_weights, nrow(source), TRUE)
  moments <- function(x, w) {
    mu <- colSums(x * (w / sum(w)))
    list(mean = mu, variance = colSums(sweep(x, 2L, mu)^2 * (w / sum(w))))
  }
  a <- moments(source, sw); b <- moments(target, tw); c <- moments(source, aw)
  sd <- sqrt((a$variance + b$variance) / 2)
  smd <- function(delta) ifelse(sd > 0, delta / sd, ifelse(delta == 0, 0, sign(delta) * Inf))
  balance <- data.frame(feature = colnames(source), before = smd(a$mean - b$mean),
    after = smd(c$mean - b$mean), reference_sd = sd, row.names = NULL)
  divergence <- brier <- NA_real_
  if (!is.null(source_prediction) || !is.null(target_prediction)) {
    if (is.null(source_prediction) || is.null(target_prediction)) stop("Both target-row prediction vectors are required for transport_div.")
    sp <- .response_p(source_prediction, nrow(target), "gaussian")
    tp <- .response_p(target_prediction, nrow(target), "gaussian")
    divergence <- sum(tw * (sp - tp)^2) / sum(tw)
  }
  if (!is.null(target_outcome)) {
    yy <- .binary(target_outcome); pp <- .prob(target_prediction, length(yy))
    if (length(yy) != nrow(target)) stop("Target outcome rows do not match.")
    brier <- sum(tw * (yy - pp)^2) / sum(tw)
  }
  list(transport_smd = balance, transport_div = divergence, target_brier = brier,
    effective_source_n = sum(aw)^2 / sum(aw^2), effective_target_n = sum(tw)^2 / sum(tw^2),
    maximum_source_weight_share = max(aw) / sum(aw),
    target = "Descriptive covariate balance and target-row prediction agreement")
}

#' Design-based decision utility, population counts and screening ratios
#' @param design A svy_learn_design.
#' @param y,p Binary outcomes and probabilities in design row order.
#' @param thresholds Thresholds in (0,1).
#' @param population_size Optional fixed population size. If supplied, counts
#'   are weighted proportions times this size; otherwise they are weight totals.
#' @param level Confidence level for pointwise intervals.
#' @param replicate_design Optional matching survey replicate design.
#' @return A long data frame. nns_svy preserves the specification's (TP+FP)/TP,
#'   interpreted as positive follow-ups per true positive. nns_population is
#'   total population / TP, appropriate when the whole population is screened.
#'   No true positives gives infinite ratios with undefined intervals. Intervals
#'   condition on predictions and are pointwise, not simultaneous bands.
#' @export
svy_decision_infer <- function(design, y, p, thresholds = seq(.01, .50, .01),
                               population_size = NULL, level = .95, replicate_design = NULL) {
  y <- .binary(y); p <- .prob(p, length(y)); .eval_level(level)
  if (length(y) != nrow(design$data)) stop("Outcome rows do not match design.")
  svy_decision(y, p, design$w, thresholds) # validates threshold schema
  if (!is.null(population_size) && (length(population_size) != 1L ||
      !is.finite(population_size) || population_size <= 0)) stop("Invalid population_size.")
  opt <- options(survey.lonely.psu = design$lonely.psu); on.exit(options(opt), add = TRUE)
  des <- .eval_design(design, replicate_design)
  out <- lapply(thresholds, function(t) {
    yes <- as.numeric(p >= t); tp <- y * yes; fp <- (1 - y) * yes
    confusion <- cbind(pop_tp = tp, pop_fp = fp, pop_tn = (1 - y) * (1 - yes), pop_fn = y * (1 - yes))
    counts <- .eval_mean(des, confusion, level, total = is.null(population_size))
    if (!is.null(population_size)) for (nm in c("estimate", "se", "ci_lower", "ci_upper"))
      counts[[nm]] <- counts[[nm]] * population_size
    benefit <- .eval_mean(des, cbind(net_benefit_svy = tp - fp * t / (1 - t),
      treat_all = y - (1 - y) * t / (1 - t), treat_none = 0), level)
    scaled <- benefit[1L, ]; scaled$metric <- "net_benefit_per_100000"
    for (nm in c("estimate", "se", "ci_lower", "ci_upper")) scaled[[nm]] <- scaled[[nm]] * 1e5
    if (sum(design$w * tp) > 0) {
      ratio <- survey::svyratio(cbind(nns_svy = yes, nns_population = rep(1, length(y))),
        matrix(tp, ncol = 1L, dimnames = list(NULL, "tp")), des)
      ratios <- .eval_table(ratio, survey::degf(des), level, "ratio; pointwise conditional interval")
      ratios$metric <- c("nns_svy", "nns_population")
    } else ratios <- data.frame(metric = c("nns_svy", "nns_population"), estimate = Inf,
      se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, df = survey::degf(des),
      method = "undefined denominator; no true positives")
    ans <- rbind(counts, benefit, scaled, ratios); ans$threshold <- t; ans
  })
  ans <- do.call(rbind, out); rownames(ans) <- NULL
  attr(ans, "population_basis") <- if (is.null(population_size))
    "Weight totals; requires population-scale weights; panel rows may represent person-waves" else
    paste("Fixed population size", population_size)
  ans
}
