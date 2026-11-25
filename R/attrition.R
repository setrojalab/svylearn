#' Specify sequential retention weighting for a panel roster
#'
#' Pass this unfitted specification as attrition in svy_fit. Retention models,
#' their preprocessing and optional wave-specific intercept stabilization are
#' re-estimated inside every inner and outer training partition. History columns
#' are lagged automatically to the preceding scheduled wave; baseline columns
#' are taken only from each subject's first row. Derived history variables must
#' themselves have been constructed without future information.
#' @param response Name of the monotone 0/1 observation indicator; baseline is 1.
#' @param history Columns to lag by one scheduled wave, including past outcomes.
#' @param baseline Optional columns taken at baseline only.
#' @param learner A single prespecified binary svy_learner; default weighted GLM.
#' @param stabilize Use a training-estimated, wave-specific intercept numerator.
#'   This preserves within-wave normalized means but changes pooled wave mass.
#' @param probability_floor Lower bound in (0,1) on fitted retention probabilities.
#'   Clipping trades bias for stability and is recorded, not a positivity remedy.
#' @param analysis Optional logical/0-1 column marking rows in the analysis target.
#'   NULL includes all observed waves. Define it before resampling.
#' @param recipe_args Arguments for training-only retention-model preprocessing.
#' @return An unfitted svy_attrition_spec. No optional engine is loaded.
#' @export
svy_attrition <- function(response, history = character(), baseline = character(),
                          learner = svy_learner(), stabilize = FALSE,
                          probability_floor = .01, analysis = NULL, recipe_args = list()) {
  if (!is.character(response) || length(response) != 1L || is.na(response) ||
      !is.character(history) || anyNA(history) || !is.character(baseline) || anyNA(baseline))
    stop("Supply response, history and baseline column names.")
  if (!inherits(learner, "svy_learner")) stop("learner must be a svy_learner specification.")
  if (!is.logical(stabilize) || length(stabilize) != 1L || is.na(stabilize)) stop("Invalid stabilize.")
  if (length(probability_floor) != 1L || !is.numeric(probability_floor) ||
      !is.finite(probability_floor) || probability_floor <= 0 || probability_floor >= 1)
    stop("probability_floor must be in (0,1).")
  if (!is.null(analysis) && (!is.character(analysis) || length(analysis) != 1L || is.na(analysis)))
    stop("analysis must name a prespecified eligibility column.")
  if (!is.list(recipe_args) || any(names(recipe_args) %in% c("data", "w"))) stop("Invalid recipe_args.")
  structure(list(response = response, history = unique(history), baseline = unique(baseline),
    learner = learner, stabilize = stabilize, probability_floor = probability_floor,
    analysis = analysis, recipe_args = recipe_args), class = "svy_attrition_spec")
}

.attrition_roster <- function(design, spec) {
  if (!inherits(design, "svy_learn_design") || !inherits(spec, "svy_attrition_spec"))
    stop("Supply a survey design and svy_attrition specification.")
  if (is.null(design$subject) || is.null(design$wave)) stop("Attrition requires subject and wave identifiers.")
  d <- design$data; n <- nrow(d); vars <- c(spec$response, spec$history, spec$baseline, spec$analysis)
  if (!all(vars %in% names(d))) stop("Unknown attrition columns.")
  time <- d[[design$wave]]
  if (!is.numeric(time) || any(!is.finite(time))) stop("Attrition wave must be finite numeric scheduled time.")
  waves <- sort(unique(time)); ids <- as.character(d[[design$subject]])
  if (length(waves) < 2L) stop("Supply baseline and at least one scheduled follow-up wave.")
  groups <- split(seq_len(n), ids)
  if (any(vapply(groups, function(i) !identical(sort(time[i]), waves), logical(1))))
    stop("Supply one row for every subject at every scheduled wave, including dropout rows.")
  response <- .binary(d[[spec$response]]); previous <- first <- integer(n)
  for (i in groups) {
    i <- i[order(time[i])]
    if (response[i[1L]] != 1 || any(diff(response[i]) > 0))
      stop("Retention must start at 1 and be monotone; intermittent missingness needs a different model.")
    if (any(design$w[i] != design$w[i[1L]])) stop("Supply constant baseline survey weights within each subject.")
    first[i] <- i[1L]; previous[i[-1L]] <- utils::head(i, -1L)
  }
  baseline <- previous == 0L
  risk <- !baseline; risk[!baseline] <- response[previous[!baseline]] == 1
  x <- data.frame(log_baseline_weight = log(design$w[first]))
  for (j in seq_along(spec$history)) x[[paste0("lag", j)]] <- d[[spec$history[j]]][pmax(1L, previous)]
  for (j in seq_along(spec$baseline)) x[[paste0("base", j)]] <- d[[spec$baseline[j]]][first]
  eligible <- response == 1
  if (!is.null(spec$analysis)) eligible <- eligible & as.logical(.binary(d[[spec$analysis]]))
  list(waves = waves, time = time, previous = previous, first = first, baseline = baseline,
       risk = risk, response = response, x = x, eligible = eligible, groups = groups)
}

#' Estimate sequential inverse-probability-of-retention weights
#'
#' The roster must retain dropout rows, have a common scheduled wave grid and
#' constant baseline survey weights within subject. The returned weights are zero
#' on unobserved rows, and baseline weights times cumulative inverse retention on
#' observed rows. Use the observed/eligible subset when constructing a new design.
#' For automatic nested fitting, pass the unfitted spec to svy_fit(attrition=spec)
#' rather than fitting this helper globally before validation.
#'
#' Identification requires sequential missing at random given the recorded history,
#' positivity and an adequate retention model. Estimated-weight uncertainty is not
#' supplied by fixed-weight survey standard errors. Agency longitudinal weights
#' should not be multiplied again by this correction. Bounded-weight sensitivity
#' does not identify an MNAR selection mechanism.
#' @param design A full-roster svy_learn_design with baseline weights.
#' @param spec An unfitted svy_attrition specification.
#' @param training Optional logical vector selecting complete training dependence
#'   groups. NULL uses all rows for a standalone fit; it is not cross-fitting.
#' @param seed Seed for optional stochastic learners; restores the caller RNG.
#' @return A svy_attrition_weights object with weights, unstabilized_weights,
#'   retention, cumulative_retention, numerator, eligible, audit and fitted model.
#' @export
#' @examples
#' d <- data.frame(id = rep(1:20, each = 2), wave = rep(0:1, 20), weight = 1,
#'                 observed = c(rbind(rep(1, 20), rep(c(1, 1, 1, 0), 5))))
#' z <- svy_design(d, "weight", subject = "id", wave = "wave")
#' a <- svy_attrition_weights(z, svy_attrition("observed"))
#' d$adjusted_weight <- a$weights
#' adjusted <- svy_design(d[a$eligible, ], "adjusted_weight", subject = "id", wave = "wave")
svy_attrition_weights <- function(design, spec, training = NULL, seed = 1L) {
  z <- .attrition_roster(design, spec); n <- nrow(design$data)
  if (is.null(training)) training <- rep(TRUE, n)
  if (!is.logical(training) || length(training) != n || anyNA(training) || !any(training))
    stop("training must select complete training groups with a logical vector.")
  if (any(vapply(split(training, .fold_cluster(design)), function(a) length(unique(a)) != 1L, logical(1))))
    stop("training must contain whole dependence groups, including complete subjects and PSUs.")
  models <- audit <- list()
  .seed(seed, {
    for (j in seq_along(z$waves)[-1L]) {
      i <- which(training & z$risk & z$time == z$waves[j])
      if (!length(i)) stop("No training subjects at risk for wave ", z$waves[j], ".")
      yy <- z$response[i]; ww <- design$w[i]
      # Deterministic no-dropout waves preserve weights exactly.
      if (all(yy == 1)) {
        f <- list(constant = 1, status = "all training subjects retained")
        rec <- NULL
      } else {
        rec <- do.call(svy_recipe, c(list(data = z$x[i, , drop = FALSE], w = ww), spec$recipe_args))
        f <- .fit_engine(bake_svy_recipe(rec, z$x[i, , drop = FALSE]), yy, ww, spec$learner)
      }
      numerator <- if (spec$stabilize) max(spec$probability_floor, sum(ww * yy) / sum(ww)) else 1
      models[[j - 1L]] <- list(recipe = rec, fit = f, numerator = numerator)
      audit[[j - 1L]] <- data.frame(wave = z$waves[j], train_at_risk = length(i),
        train_retained = sum(yy), numerator = numerator, model_status = f$status)
    }
  })
  model <- structure(list(spec = spec, waves = z$waves, models = models,
    training_clusters = unique(.fold_cluster(design)[training]), audit = do.call(rbind, audit)),
    class = "svy_attrition_model")
  out <- .attrition_apply(model, design, z); out$model <- model; out
}

.attrition_apply <- function(model, design, z = NULL) {
  if (is.null(z)) z <- .attrition_roster(design, model$spec)
  if (!identical(z$waves, model$waves)) stop("Prediction roster must use the fitted scheduled wave grid.")
  n <- nrow(design$data); p <- raw <- numerator <- rep(NA_real_, n)
  p[z$baseline] <- raw[z$baseline] <- numerator[z$baseline] <- 1
  audit <- model$audit
  for (j in seq_along(z$waves)[-1L]) {
    i <- which(z$risk & z$time == z$waves[j]); f <- model$models[[j - 1L]]
    if (length(i)) {
      raw[i] <- if (is.null(f$recipe)) f$fit$constant else
        .predict_engine(f$fit, bake_svy_recipe(f$recipe, z$x[i, , drop = FALSE]))
      p[i] <- pmax(model$spec$probability_floor, raw[i]); numerator[i] <- f$numerator
    }
    audit$applied_at_risk[j - 1L] <- length(i)
    audit$clipped[j - 1L] <- sum(raw[i] < model$spec$probability_floor)
  }
  cumulative <- cumulative_num <- rep(NA_real_, n)
  cumulative[z$baseline] <- cumulative_num[z$baseline] <- 1
  for (j in seq_along(z$waves)[-1L]) {
    i <- which(z$risk & z$time == z$waves[j]); previous <- z$previous[i]
    cumulative[i] <- cumulative[previous] * p[i]
    cumulative_num[i] <- cumulative_num[previous] * numerator[i]
  }
  weights <- unstabilized <- numeric(n); observed <- z$response == 1
  unstabilized[observed] <- design$w[observed] / cumulative[observed]
  weights[observed] <- unstabilized[observed] * cumulative_num[observed]
  if (any(!is.finite(weights[observed])) || any(weights[observed] <= 0))
    stop("Cumulative retention weights overflowed or underflowed; inspect positivity and wave count.")
  structure(list(weights = weights, unstabilized_weights = unstabilized, retention = p,
    raw_retention = raw, cumulative_retention = cumulative, numerator = numerator,
    eligible = z$eligible, audit = audit,
    interpretation = "Estimated sequential-MAR retention weights; no MNAR identification or fixed-weight nuisance variance"),
    class = "svy_attrition_weights")
}

.attrition_observed <- function(d, weights, eligible) {
  i <- which(eligible)
  if (!length(i)) stop("No observed eligible outcomes in a training or evaluation partition.")
  out <- .subset_design(d, i); out$w <- weights[i]; out$data[[out$weights]] <- out$w
  # Agency baseline replicate weights do not include fitted retention uncertainty.
  out$repweights <- NULL; out$replicate_df <- NULL
  out
}
