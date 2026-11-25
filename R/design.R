#' Specify a probability sample and optional panel nesting
#'
#' Stores design metadata without loading the survey package. Nested subject
#' identifiers must be globally unique. Repeated subject waves must remain in
#' one highest-stage PSU. Supply the appropriate wave or attrition weights;
#' their validity is an assumption, not something inferred from the data.
#' @param data A data frame.
#' @param weights Name or one-sided formula for the sampling weight.
#' @param ids Names or formula for ordered sampling stages; NULL for independent units.
#' @param strata Names or formula for strata, optionally at each stage.
#' @param fpc Names or formula for finite population corrections.
#' @param subject,wave,household Optional column names for panel identifiers.
#' @param repweights Optional numeric replicate-weight matrix.
#' @param type Replicate type understood by survey::svrepdesign.
#' @param scale,rscales,mse,combined.weights Replicate variance conventions.
#' @param lonely.psu Singleton policy: fail, adjust, average, certainty or remove.
#' @param fold_groups Optional column names defining additional dependence groups
#'   (sites, spatial blocks, graph components). Groups and sampling PSUs are joined
#'   transitively for validation only; inference retains the supplied sampling design.
#' @param replicate_df Optional agency-specified positive degrees of freedom for
#'   supplied replicate weights. NULL uses survey's rank-based default.
#' @return An object of class svy_learn_design.
#' @export
#' @examples
#' d <- svy_example(200)
#' design <- svy_design(d, "weight", "psu", "stratum")
#' svy_diagnostics(design)
svy_design <- function(data, weights, ids = NULL, strata = NULL, fpc = NULL,
                       subject = NULL, wave = NULL, household = NULL,
                       repweights = NULL, type = "other", scale = NULL,
                       rscales = NULL, mse = TRUE, combined.weights = TRUE,
                       lonely.psu = "fail", fold_groups = NULL, replicate_df = NULL) {
  if (!is.data.frame(data) || !nrow(data)) stop("data must be a nonempty data frame.")
  weights <- .cols(weights); ids <- .cols(ids); strata <- .cols(strata); fpc <- .cols(fpc)
  fold_groups <- .cols(fold_groups)
  vars <- unique(c(weights, ids, strata, fpc, subject, wave, household, fold_groups))
  if (length(weights) != 1L || !all(vars %in% names(data))) stop("Unknown or invalid design columns.")
  if (anyNA(data[vars])) stop("Design variables cannot contain missing values.")
  w <- .positive(data[[weights]], nrow(data))
  skey <- .key(data, utils::head(strata, 1L))
  cluster <- if (length(ids)) .key(data, c(utils::head(strata, 1L), ids[1L])) else as.character(seq_len(nrow(data)))
  if (!length(ids) && !is.null(subject)) cluster <- .key(data, subject)
  if (any(vapply(split(skey, cluster), function(z) length(unique(z)) > 1L, logical(1))))
    stop("A stable sampling cluster cannot cross strata.")
  for (v in c(subject, household)) {
    nesting <- split(cluster, data[[v]])
    if (any(vapply(nesting, function(z) length(unique(z)) > 1L, logical(1))))
      stop(v, " crosses highest-stage PSUs; specify a stable common cluster before splitting.")
  }
  if (!is.null(subject) && !is.null(wave) && anyDuplicated(data[c(subject, wave)]))
    stop("Duplicate subject-wave records.")
  lonely.psu <- match.arg(lonely.psu, c("fail", "adjust", "average", "certainty", "remove"))
  if (!is.null(replicate_df) && (is.null(repweights) || length(replicate_df) != 1L ||
      !is.numeric(replicate_df) || !is.finite(replicate_df) || replicate_df <= 0 || replicate_df > ncol(repweights)))
    stop("replicate_df requires supplied replicate weights and a positive finite value.")
  if (!is.null(repweights)) {
    if (!is.matrix(repweights) || !is.numeric(repweights) || nrow(repweights) != nrow(data) ||
        ncol(repweights) < 2L || any(!is.finite(repweights)) || any(repweights < 0) ||
        any(colSums(repweights) <= 0)) stop("Invalid replicate weights.")
    if (type == "other" && (is.null(scale) || is.null(rscales)))
      stop("Custom replicate weights require explicit scale and rscales.")
  }
  structure(list(data = data, weights = weights, w = w, ids = ids, strata = strata,
    fpc = fpc, subject = subject, wave = wave, household = household,
    cluster = cluster, stratum = skey, repweights = repweights, type = type,
    scale = scale, rscales = rscales, mse = mse, combined.weights = combined.weights,
    lonely.psu = lonely.psu, fold_groups = fold_groups, replicate_df = replicate_df,
    fold_cluster = if (length(fold_groups)) .dependency_groups(cluster, data[fold_groups]) else NULL),
    class = "svy_learn_design")
}

#' Survey weight and cluster diagnostics
#' @param design A svy_learn_design.
#' @param y Optional binary outcome.
#' @return A named list of sample, weight and cluster summaries.
#' @export
svy_diagnostics <- function(design, y = NULL) {
  w <- design$w
  ct <- table(design$stratum[!duplicated(design$cluster)])
  out <- list(n = length(w), psus = length(unique(design$cluster)), strata = length(ct),
    singleton_strata = sum(ct == 1L), design_df = sum(ct - 1L),
    weight_sum = sum(w), effective_n = sum(w)^2 / sum(w^2),
    weight_cv = sqrt(mean((w / mean(w) - 1)^2)), weight_max_ratio = max(w) / min(w))
  if (!is.null(y)) {
    y <- .binary(y); if (length(y) != length(w)) stop("Outcome length mismatch.")
    ew <- w[y == 1]
    out$events <- sum(y); out$prevalence <- sum(w * y) / sum(w)
    out$effective_events <- if (length(ew)) sum(ew)^2 / sum(ew^2) else 0
  }
  out
}

#' Convert to a survey design on demand
#' @param design A svy_learn_design.
#' @return A survey.design or svyrep.design object.
#' @export
as_survey_design <- function(design) {
  .need("survey")
  form <- function(x, default = "~1") if (length(x)) stats::reformulate(x) else stats::as.formula(default)
  if (is.null(design$repweights)) {
    ids <- design$ids
    if (!length(ids) && !is.null(design$subject)) ids <- design$subject
    survey::svydesign(ids = form(ids), strata = if (length(design$strata)) form(design$strata) else NULL,
      weights = form(design$weights), fpc = if (length(design$fpc)) form(design$fpc) else NULL,
      data = design$data, nest = TRUE)
  } else {
    rd <- survey::svrepdesign(data = design$data, weights = form(design$weights),
      repweights = design$repweights, type = design$type, scale = design$scale,
      rscales = design$rscales, mse = design$mse, combined.weights = design$combined.weights,
      degf = design$replicate_df)
    attr(rd, "svylearn_replicate_df") <- design$replicate_df
    rd
  }
}

.survey_domain <- function(design, domain) {
  specified <- attr(design, "svylearn_replicate_df")
  out <- design[domain, ]
  # survey recomputes rank-based degrees of freedom on subsetting. Retain an
  # explicitly supplied agency convention for parent-design domain inference.
  if (!is.null(specified) && inherits(out, "svyrep.design")) {
    out$degf <- specified; attr(out, "svylearn_replicate_df") <- specified
  }
  out
}

#' Construct survey replicate weights
#' @param design A svy_learn_design.
#' @param type Replication scheme, such as bootstrap, BRR, Fay, JK1 or JKn.
#' @param replicates Number of bootstrap replicates.
#' @param seed Random seed, restored on exit.
#' @param ... Further arguments to survey::as.svrepdesign.
#' @return A survey replicate design.
#' @export
svy_replicates <- function(design, type = "bootstrap", replicates = 500L, seed = 1L, ...) {
  .need("survey")
  args <- c(list(design = as_survey_design(design), type = type), list(...))
  if (type %in% c("bootstrap", "subbootstrap", "mrbbootstrap")) args$replicates <- replicates
  .seed(seed, do.call(survey::as.svrepdesign, args))
}

.subset_design <- function(d, i) {
  d$data <- d$data[i, , drop = FALSE]; d$w <- d$w[i]
  d$cluster <- d$cluster[i]; d$stratum <- d$stratum[i]
  if (!is.null(d$fold_cluster)) d$fold_cluster <- d$fold_cluster[i]
  if (!is.null(d$repweights)) d$repweights <- d$repweights[i, , drop = FALSE]
  d
}

#' Allocate whole sampling clusters to validation folds
#' @param design A svy_learn_design.
#' @param v,repeats Number of folds and repeats.
#' @param seed Random seed.
#' @return An integer matrix, one column per repeat. Assignment uses design
#'   variables only, never outcomes. With fewer than v PSUs in a stratum that
#'   stratum cannot be represented in every validation fold.
#' @export
svy_split <- function(design, v = 5L, repeats = 1L, seed = 1L) {
  if (v < 2L || v != as.integer(v) || repeats < 1L || repeats != as.integer(repeats))
    stop("v and repeats must be positive integers; v must be at least two.")
  cluster <- .fold_cluster(design)
  keys <- unique(cluster); first <- match(keys, cluster)
  if (length(keys) < v) stop("Fewer independent clusters than requested folds.")
  crosses <- any(vapply(split(design$stratum, cluster), function(x) length(unique(x)) > 1L, logical(1)))
  groups <- if (crosses) list(seq_along(keys)) else split(seq_along(keys), design$stratum[first])
  idx <- match(cluster, keys)
  mass <- as.numeric(rowsum(design$w, idx, reorder = FALSE))
  .seed(seed, {
    out <- matrix(0L, length(idx), repeats)
    for (r in seq_len(repeats)) {
      allocation <- integer(length(keys)); total <- numeric(v)
      for (gg in sample.int(length(groups))) {
        g <- groups[[gg]]; g <- g[sample.int(length(g))]
        local <- integer(v)
        for (j in g) {
          candidates <- which(local == min(local))
          f <- candidates[which.min(total[candidates])]
          allocation[j] <- f; local[f] <- local[f] + 1L; total[f] <- total[f] + mass[j]
        }
      }
      out[, r] <- allocation[idx]
    }
    colnames(out) <- paste0("repeat", seq_len(repeats)); out
  })
}

.fold_cluster <- function(d) if (is.null(d$fold_cluster)) d$cluster else d$fold_cluster
.dependency_groups <- function(cluster, groups) {
  keys <- unique(cluster); index <- match(cluster, keys); parent <- seq_along(keys)
  root <- function(i) {
    while (parent[i] != i) { parent[i] <<- parent[parent[i]]; i <- parent[i] }
    i
  }
  for (v in groups) for (g in split(index, as.character(v))) {
    g <- unique(g)
    if (length(g) > 1L) {
      a <- root(g[1L]); for (b in g[-1L]) parent[root(b)] <- a
    }
  }
  vapply(index, root, integer(1))
}
