#' Refit the complete nested workflow for survey replicate weights
#'
#' Repeats training, tuning and calibration for each replicate-weight column.
#' Unlike fixed-score inference, this captures variation in the complete
#' statistic under the supplied replication scheme. Its inferential validity
#' still depends on the sampling design, replicate construction and learning
#' procedure. Failed replicates are retained in the audit and prevent an SE.
#' @param fit A svy_learn_fit.
#' @param replicate_design A survey replicate design with the same row ordering.
#'   For attrition fits this must describe the full baseline roster, including
#'   dropout rows, with replicate weights constant within subjects. Retention
#'   models and their preprocessing are refitted in every nested partition.
#' @return A list with replicate metrics, statuses and the variance matrix when
#'   every replicate succeeds. Computation is sequential and deterministic.
#' @export
svy_refit_replicates <- function(fit, replicate_design) {
  .need("survey")
  ww <- as.matrix(stats::weights(replicate_design, type = "analysis"))
  source <- if (inherits(fit, "svy_attrition_fit")) fit$attrition$roster else fit$design
  if (nrow(ww) != nrow(source$data)) stop("Replicate row count mismatch; attrition fits require the full baseline roster.")
  metric_names <- if (inherits(fit, "svy_response_fit") || inherits(fit, "svy_attrition_fit")) colnames(fit$metrics) else
    c("log_loss", "brier", "auc", "average_precision")
  estimates <- matrix(NA_real_, ncol(ww), length(metric_names), dimnames = list(NULL, metric_names))
  status <- rep("ok", ncol(ww))
  for (r in seq_len(ncol(ww))) {
    tryCatch({
      active <- which(ww[, r] > 0)
      d <- .subset_design(source, active)
      d$w <- ww[active, r]; d$data[[d$weights]] <- d$w; d$repweights <- NULL
      conf <- fit$configuration
      ff <- do.call(svy_fit, c(list(design = d), conf))
      estimates[r, ] <- colMeans(ff$metrics)[metric_names]
    }, error = function(e) status[r] <<- conditionMessage(e))
  }
  variance <- NULL
  if (all(status == "ok") && all(is.finite(estimates))) {
    variance <- survey::svrVar(estimates, scale = replicate_design$scale,
      rscales = replicate_design$rscales, mse = replicate_design$mse,
      coef = colMeans(fit$metrics)[metric_names])
  }
  list(estimates = estimates, status = status, variance = variance,
    target = "Replicated complete nested CV statistic under the supplied replicate scheme")
}
