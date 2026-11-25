# Shared inference helpers. Optional survey is loaded only at evaluation time.
.eval_level <- function(level) {
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1)
    stop("level must be strictly between zero and one.")
  level
}
.eval_design <- function(design, replicate_design = NULL) {
  .need("survey")
  if (!inherits(design, "svy_learn_design")) stop("Supply a svy_learn_design.")
  if (is.null(replicate_design)) return(as_survey_design(design))
  rd <- replicate_design
  if (!inherits(rd, "svyrep.design") || nrow(rd$variables) != nrow(design$data))
    stop("replicate_design must be a survey replicate design in the same row order.")
  if (!isTRUE(all.equal(as.numeric(stats::weights(rd, type = "sampling")),
                        design$w, tolerance = 1e-10)))
    stop("Replicate design sampling weights do not match design weights.")
  keys <- unique(c(design$ids, design$strata, design$subject, design$wave))
  if (length(keys) && (!all(keys %in% names(rd$variables)) ||
      !all(vapply(keys, function(k) identical(as.character(rd$variables[[k]]),
        as.character(design$data[[k]])), logical(1)))))
    stop("Replicate design identifiers or row order do not match.")
  rd
}
.eval_table <- function(est, df, level, method) {
  .eval_level(level)
  val <- stats::coef(est); se <- sqrt(pmax(0, diag(stats::vcov(est))))
  labels <- names(val)
  if (is.null(labels)) labels <- if (!is.null(dim(val)) && ncol(val) == 1L) rownames(val) else NULL
  if (is.null(labels)) labels <- paste0("value", seq_along(val))
  q <- if (is.finite(df) && df > 0) stats::qt((1 + level) / 2, df) else NA_real_
  data.frame(metric = labels, estimate = as.numeric(val), se = as.numeric(se),
    ci_lower = as.numeric(val - q * se), ci_upper = as.numeric(val + q * se),
    df = df, method = method, row.names = NULL)
}
.eval_mean <- function(des, values, level = .95, total = FALSE) {
  if (is.null(dim(values))) values <- matrix(values, ncol = 1L)
  if (!is.numeric(values) || nrow(values) != nrow(des$variables) ||
      any(!is.finite(values))) stop("Evaluation values must be finite and match design rows.")
  if (is.null(colnames(values))) colnames(values) <- paste0("value", seq_len(ncol(values)))
  est <- if (total) survey::svytotal(values, des) else survey::svymean(values, des)
  .eval_table(est, survey::degf(des), level,
    if (inherits(des, "svyrep.design")) "replicate; fixed predictions" else "Taylor; fixed predictions")
}
.eval_replicate <- function(des, fun, level = .95) {
  if (!inherits(des, "svyrep.design")) stop("Nonlinear replication requires a survey replicate design.")
  est <- survey::withReplicates(des, function(w, data) {
    ans <- fun(w)
    if (any(!is.finite(ans))) stop("A replicate statistic is undefined; no replicates were silently dropped.")
    ans
  })
  .eval_table(est, survey::degf(des), level, "replicate; statistic recomputed")
}
.eval_losses <- function(y, p, family) {
  if (family == "binomial") {
    y <- .binary(y); p <- .prob(p, length(y)); pc <- .clip(p)
    return(cbind(log_loss = -(y * log(pc) + (1 - y) * log1p(-pc)), brier = (y - p)^2))
  }
  y <- .response_y(y, family); p <- .response_p(p, length(y), family, levels(y))
  .response_losses(y, p, family)
}
.eval_calibration <- function(des, y, p, family, level) {
  if (!family %in% c("binomial", "gaussian")) return(NULL)
  if (length(unique(y)) < 2L || length(unique(p)) < 2L) return(NULL)
  des$variables$.svye_y <- y
  des$variables$.svye_p <- if (family == "binomial") stats::qlogis(.clip(p)) else p
  if (family == "binomial") {
    joint <- survey::svyglm(.svye_y ~ .svye_p, des, family = stats::quasibinomial())
    offset <- survey::svyglm(.svye_y ~ 1 + offset(.svye_p), des, family = stats::quasibinomial())
  } else {
    joint <- survey::svyglm(.svye_y ~ .svye_p, des, family = stats::gaussian())
    offset <- survey::svyglm(.svye_y ~ 1 + offset(.svye_p), des, family = stats::gaussian())
  }
  if (!isTRUE(joint$converged) || !isTRUE(offset$converged) ||
      any(!is.finite(stats::coef(joint))) || any(!is.finite(stats::coef(offset)))) return(NULL)
  jt <- .eval_table(joint, joint$df.residual, level, "joint calibration regression")
  ot <- .eval_table(offset, offset$df.residual, level, "offset calibration regression")
  jt$metric <- c("cal_joint_intercept", "cal_slope"); ot$metric <- "cal_intercept"
  rbind(ot, jt)
}
