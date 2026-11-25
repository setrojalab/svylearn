#' Generate a reproducible survey or panel teaching dataset
#' @param n Approximate cross-sectional sample size, at least 100.
#' @param waves Number of waves; one gives a cross-sectional sample.
#' @param seed Random seed.
#' @return A synthetic data frame. It is not NHANES, HRS or PSID data.
#' @export
svy_example <- function(n = 1000L, waves = 1L, seed = 17L) {
  if (n < 100 || n != as.integer(n) || waves < 1 || waves != as.integer(waves)) stop("Invalid n or waves.")
  .seed(seed, {
    stratum <- rep(seq_len(10), length.out = n)
    psu <- rep(rep(seq_len(4), each = 10), length.out = n)
    cluster <- interaction(stratum, psu, drop = TRUE)
    x1 <- stats::rnorm(n); x2 <- stats::rnorm(n)
    re <- stats::rnorm(nlevels(cluster), sd = .4)[as.integer(cluster)]
    weight <- exp(.5 * x1)
    base <- data.frame(subject = seq_len(n), stratum = stratum, psu = psu,
      weight = weight, x1 = x1, x2 = x2, re = re)
    d <- base[rep(seq_len(n), each = waves), ]; d$wave <- rep(seq_len(waves), n)
    d$y <- stats::rbinom(nrow(d), 1, stats::plogis(-2 + .7 * d$x1 - .5 * d$x2 + d$re + .15 * (d$wave - 1)))
    d$re <- NULL; rownames(d) <- NULL; d
  })
}

#' Simulate survey sampling around a finite population
#'
#' This experiment evaluates fixed-score inference with stratified simple
#' random sampling of PSUs and a census of units in selected PSUs. The supplied
#' prediction vector is fixed before sampling. It does not validate the coverage
#' of a learned model or of nonresponse-adjusted estimators.
#' @param population Data frame containing outcome, score, PSU and stratum columns.
#' @param outcome,prediction,psu,stratum Column names.
#' @param sampled_psus Number of PSUs sampled per stratum, at least two.
#' @param replicates Number of repeated probability samples.
#' @param seed Random seed.
#' @return A list with replicate estimates, truth and bias, RMSE, coverage and MCSE.
#' @export
svy_plasmode <- function(population, outcome, prediction, psu, stratum,
                         sampled_psus = 2L, replicates = 500L, seed = 1L) {
  .need("survey")
  y <- .binary(population[[outcome]]); p <- .prob(population[[prediction]], length(y))
  if (sampled_psus < 2 || sampled_psus != as.integer(sampled_psus) || replicates < 2) stop("Invalid replication settings.")
  keys <- .key(population, c(stratum, psu)); first <- !duplicated(keys)
  groups <- split(keys[first], population[[stratum]][first]); counts <- lengths(groups)
  if (any(counts <= sampled_psus)) stop("Each stratum must contain more PSUs than sampled_psus.")
  truth <- svy_metrics(y, p)[c("log_loss", "brier")]
  .seed(seed, {
    out <- vector("list", replicates)
    for (r in seq_len(replicates)) {
      selected <- unlist(lapply(groups, function(z) z[sample.int(length(z), sampled_psus)]), use.names = FALSE)
      take <- keys %in% selected; d <- population[take, , drop = FALSE]
      d$.svyl_fpc <- as.numeric(counts[as.character(d[[stratum]])])
      d$.svyl_w <- d$.svyl_fpc / sampled_psus
      des <- svy_design(d, ".svyl_w", psu, stratum, fpc = ".svyl_fpc")
      z <- svy_infer(des, d[[outcome]], d[[prediction]])$scores
      z$replicate <- r; z$truth <- as.numeric(truth[z$metric]); out[[r]] <- z
    }
    records <- do.call(rbind, out)
    summary <- do.call(rbind, lapply(split(records, records$metric), function(z) {
      b <- z$estimate - z$truth; cvr <- mean(z$lower <= z$truth & z$upper >= z$truth)
      data.frame(metric = z$metric[1], replicates = nrow(z), bias = mean(b),
        mcse_bias = stats::sd(b) / sqrt(nrow(z)), rmse = sqrt(mean(b^2)),
        coverage = cvr, mcse_coverage = sqrt(cvr * (1 - cvr) / nrow(z)))
    }))
    list(replicates = records, truth = truth, summary = summary)
  })
}
