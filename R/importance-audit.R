.permutation_rows <- function(d) {
  # Panel predictors travel as whole subject trajectories among compatible schedules.
  if (is.null(d$subject) || is.null(d$wave)) return(split(seq_len(nrow(d$data)), d$cluster))
  subject <- split(seq_len(nrow(d$data)), as.character(d$data[[d$subject]]))
  subject <- lapply(subject, function(i) i[order(d$data[[d$wave]][i])])
  schedule <- vapply(subject, function(i) .key(data.frame(
    cluster = d$cluster[i[1L]], schedule = paste(.key(d$data[i, , drop = FALSE], d$wave), collapse = "|")),
    c("cluster", "schedule")), character(1))
  lapply(split(seq_along(subject), schedule), function(i) do.call(cbind, subject[i]))
}
.permute_feature <- function(x, groups) {
  out <- x
  for (g in groups) {
    if (is.matrix(g)) out[as.vector(g)] <- x[as.vector(g[, sample.int(ncol(g)), drop = FALSE])] else
      out[g] <- x[g[sample.int(length(g))]]
  }
  out
}

#' Conditional design intervals for within-PSU permutation importance
#' @param fit A fitted svy_learn_fit with retained outer models.
#' @param repeats Number of Monte Carlo permutations per CV repeat.
#' @param seed Random seed, restored on exit.
#' @param level Confidence level.
#' @param replicate_design Optional matching survey replicate design, including
#'   a cluster-bootstrap design when bootstrap intervals are desired.
#' @return vip_psu data frame and per-draw audit. The interval treats fitted models,
#'   realized permutations and their row-level loss differences as fixed. The
#'   separate Monte Carlo SE does not represent training uncertainty. With subject
#'   and wave metadata, only whole trajectories with identical wave schedules
#'   are permuted within the same PSU. Unmatched trajectories remain unchanged.
#' @export
svy_importance_infer <- function(fit, repeats = 10L, seed = 1L, level = .95,
                                 replicate_design = NULL) {
  if (!inherits(fit, "svy_learn_fit")) stop("Supply a fitted svy_learn_fit.")
  if (length(repeats) != 1L || !is.finite(repeats) || repeats < 1L || repeats != as.integer(repeats))
    stop("repeats must be a positive integer.")
  .eval_level(level)
  d <- fit$design; y <- d$data[[fit$outcome]]; generalized <- inherits(fit, "svy_response_fit")
  family <- if (generalized) fit$family else "binomial"
  opt <- options(survey.lonely.psu = d$lonely.psu); on.exit(options(opt), add = TRUE)
  des <- .eval_design(d, replicate_design); groups <- .permutation_rows(d)
  .seed(seed, {
    n <- nrow(d$data); nr <- ncol(fit$folds); v <- max(fit$folds)
    delta <- matrix(0, n, length(fit$predictors), dimnames = list(NULL, fit$predictors))
    audit <- list()
    for (r in seq_len(nr)) {
      base <- if (family == "multinomial") matrix(fit$predictions[, , r], n,
        length(fit$levels), dimnames = list(NULL, fit$levels)) else fit$predictions[, r]
      base_loss <- .eval_losses(y, base, family)[, 1L]
      for (j in seq_len(repeats)) for (variable in fit$predictors) {
        nd <- d$data; nd[[variable]] <- .permute_feature(nd[[variable]], groups)
        pred <- if (family == "multinomial") matrix(0, n, length(fit$levels),
          dimnames = list(NULL, fit$levels)) else numeric(n)
        for (k in seq_len(v)) {
          i <- which(fit$folds[, r] == k); f <- fit$outer_fits[[(r - 1L) * v + k]]
          pp <- if (generalized) .response_pipeline_predict(f, nd[i, , drop = FALSE]) else
            .pipeline_predict(f, nd[i, , drop = FALSE])
          if (family == "multinomial") pred[i, ] <- pp else pred[i] <- pp
        }
        dd <- .eval_losses(y, pred, family)[, 1L] - base_loss
        delta[, variable] <- delta[, variable] + dd / (repeats * nr)
        audit[[length(audit) + 1L]] <- data.frame(feature = variable, repeat_id = r,
          permutation = j, loss_increase = sum(d$w * dd) / sum(d$w))
      }
    }
    audit <- do.call(rbind, audit)
    table <- .eval_mean(des, delta, level); names(table)[1L] <- "feature"
    table$vip_psu <- table$estimate
    table$monte_carlo_se <- vapply(table$feature, function(variable) {
      a <- audit[audit$feature == variable, ]
      draws <- tapply(a$loss_increase, a$permutation, mean)
      if (length(draws) > 1L) stats::sd(draws) / sqrt(length(draws)) else NA_real_
    }, numeric(1))
    movable <- sum(vapply(groups, function(g) if (is.matrix(g))
      if (ncol(g) > 1L) length(g) else 0L else if (length(g) > 1L) length(g) else 0L, integer(1)))
    list(vip_psu = table, draws = audit, permutable_row_fraction = movable / n,
      target = "Conditional associational importance; fixed trained models and realized permutations",
      permutation_unit = if (is.null(d$subject) || is.null(d$wave)) "rows within PSU" else
        "whole subjects with compatible wave schedules within PSU")
  })
}
