attrition_fixture <- function(n = 240L) {
  set.seed(883)
  x <- rnorm(n); y0 <- x + rnorm(n); y1 <- .7 * y0 + .3 * x + rnorm(n)
  y2 <- .7 * y1 + .3 * x + rnorm(n)
  r1 <- rbinom(n, 1, plogis(1 - .5 * y0 + .2 * x))
  r2 <- r1 * rbinom(n, 1, plogis(1 - .5 * y1 + .2 * x))
  d <- data.frame(id = rep(seq_len(n), each = 3), wave = rep(0:2, n),
    psu = rep(rep(seq_len(n / 8), each = 8), each = 3), w = rep(exp(.2 * x), each = 3),
    x = rep(x, each = 3), y = as.vector(t(cbind(y0, y1, y2))),
    observed = as.vector(t(cbind(1, r1, r2))))
  d$y[d$observed == 0] <- NA_real_; d$target <- d$wave == 2
  svy_design(d, "w", "psu", subject = "id", wave = "wave")
}

test_that("sequential weights match independent weighted logistic calculations", {
  z <- attrition_fixture(); d <- z$data
  spec <- svy_attrition("observed", history = c("y", "x"), analysis = "target")
  a <- svy_attrition_weights(z, spec)
  p <- rep(NA_real_, nrow(d)); p[d$wave == 0] <- 1
  for (t in 1:2) {
    i <- which(d$wave == t & d$observed[pmax(1, seq_len(nrow(d)) - 1)] == 1)
    oracle <- data.frame(r = d$observed[i], lag_y = d$y[i - 1], lag_x = d$x[i - 1], logw = log(d$w[i]))
    fit <- stats::glm(r ~ lag_y + lag_x + logw, data = oracle,
      weights = d$w[i] / mean(d$w[i]), family = quasibinomial())
    p[i] <- pmax(.01, as.numeric(stats::predict(fit, type = "response")))
  }
  expect_equal(a$retention, p, tolerance = 1e-8)
  keep <- which(d$observed == 1 & d$wave == 2)
  expect_equal(a$weights[keep], d$w[keep] / (p[keep - 1] * p[keep]), tolerance = 1e-8)
  expect_equal(a$weights[d$observed == 0], rep(0, sum(d$observed == 0)))
  expect_identical(a$eligible, d$observed == 1 & d$target)
  expect_equal(a$weights[d$wave == 0], d$w[d$wave == 0])
})

test_that("stabilization uses training risk sets and has an explicit wave target", {
  z <- attrition_fixture(); d <- z$data
  a <- svy_attrition_weights(z, svy_attrition("observed", "y"))
  b <- svy_attrition_weights(z, svy_attrition("observed", "y", stabilize = TRUE))
  keep <- which(d$observed == 1 & d$wave == 2)
  ratio <- b$weights[keep] / a$weights[keep]
  expect_equal(ratio, rep(prod(b$model$audit$numerator), length(keep)))
  expect_equal(weighted.mean(d$y[keep], a$weights[keep]), weighted.mean(d$y[keep], b$weights[keep]))
  expect_equal(b$unstabilized_weights, a$weights)
  expect_lt(prod(b$model$audit$numerator), 1)
  clipped <- svy_attrition_weights(z, svy_attrition("observed", "y", probability_floor = .95))
  expect_true(any(clipped$audit$clipped > 0))
  expect_true(all(clipped$retention[is.finite(clipped$retention)] >= .95))
})

test_that("complete rosters and monotone retention are required", {
  z <- attrition_fixture(); spec <- svy_attrition("observed", "y")
  expect_error(svy_attrition_weights(svylearn:::.subset_design(z, 2:nrow(z$data)), spec), "every scheduled wave")
  bad <- z; bad$data$observed[1:3] <- c(1, 0, 1)
  expect_error(svy_attrition_weights(bad, spec), "monotone")
  bad <- z; bad$w[2] <- 10
  expect_error(svy_attrition_weights(bad, spec), "constant baseline")
  tr <- rep(TRUE, nrow(z$data)); tr[1] <- FALSE
  expect_error(svy_attrition_weights(z, spec, tr), "whole dependence groups")
  expect_error(svy_attrition("observed", probability_floor = 0), "in \\(0,1\\)")
  expect_error(svy_attrition("observed", stabilize = NA), "Invalid stabilize")
  expect_error(svy_attrition_weights(z, svy_attrition("unknown")), "Unknown attrition")
  no_panel <- z; no_panel$subject <- NULL
  expect_error(svy_attrition_weights(no_panel, spec), "subject and wave")
})

test_that("no attrition preserves baseline weights without smoothing bias", {
  z <- attrition_fixture(); z$data$observed <- 1
  a <- svy_attrition_weights(z, svy_attrition("observed"))
  expect_identical(a$weights, z$w)
  expect_equal(a$retention, rep(1, nrow(z$data)))
  expect_true(all(a$audit$model_status == "all training subjects retained"))
})

test_that("current and future values cannot enter an earlier retention model", {
  z <- attrition_fixture(); spec <- svy_attrition("observed", history = "y", baseline = "x")
  a <- svy_attrition_weights(z, spec)
  later <- z; later$data$y[later$data$wave == 2] <- 1e6
  later$data$x[later$data$wave > 0] <- -1e6
  b <- svy_attrition_weights(later, spec)
  expect_identical(a$retention, b$retention)
  tr <- svy_split(z, 3, seed = 9)[, 1] != 1
  changed <- z; changed$data$observed[!tr] <- 1; changed$data$y[!tr] <- 200
  a <- svy_attrition_weights(z, spec, tr)
  b <- svy_attrition_weights(changed, spec, tr)
  expect_identical(a$model, b$model)
})

test_that("nested attrition training excludes held-out PSUs at both levels", {
  z <- attrition_fixture(); spec <- svy_attrition("observed", c("y", "x"), analysis = "target")
  fit <- svy_fit(z, "y", "x", family = "gaussian", v = 3, inner_v = 2, attrition = spec, seed = 18)
  te <- fit$attrition$roster_folds[, 1] == 1
  changed <- z; changed$data$y[te & changed$data$observed == 1] <- 100
  other <- svy_fit(changed, "y", "x", family = "gaussian", v = 3, inner_v = 2, attrition = spec, seed = 18)
  expect_identical(fit$outer_fits[[1]], other$outer_fits[[1]])
  expect_equal(fit$predictions[fit$folds[, 1] == 1, ], other$predictions[other$folds[, 1] == 1, ])
  for (k in 1:3) {
    tr <- which(fit$attrition$roster_folds[, 1] != k); f <- fit$outer_fits[[k]]
    expect_false(any(z$cluster[!seq_len(nrow(z$data)) %in% tr] %in% f$attrition_model$training_clusters))
    for (j in 1:2) expect_setequal(f$inner_attrition_audit[[j]]$training_clusters,
      unique(z$cluster[tr][f$roster_inner_folds != j]))
  }
  expect_true(all(is.finite(fit$attrition$weights)))
  expect_equal(fit$metrics[1, ], svy_metrics(fit$design$data$y, fit$predictions[, 1], fit$attrition$weights[, 1], "gaussian"))
  path <- tempfile(); saveRDS(fit, path)
  expect_identical(predict(readRDS(path), fit$design$data), predict(fit, fit$design$data))
  skip_if_not_installed("survey")
  expect_match(svy_eval(fit)$prediction_source, "retention", fixed = TRUE)
})

test_that("nested attrition supports all four response schemas", {
  z <- attrition_fixture(); observed <- z$data$observed == 1
  spec <- svy_attrition("observed", "x", analysis = "target")
  for (family in c("binomial", "poisson", "multinomial")) {
    d <- z
    d$data$outcome <- switch(family, binomial = as.numeric(d$data$y > 0),
      poisson = round(exp(pmin(2, d$data$y / 2))),
      multinomial = cut(d$data$y, c(-Inf, -.5, .5, Inf), labels = c("a", "b", "c")))
    f <- svy_fit(d, "outcome", "x", family = family, calibration = "none", v = 3, inner_v = 2, attrition = spec)
    expect_true(all(is.finite(f$predictions)))
    expect_equal(nrow(f$design$data), sum(observed & z$data$target))
    expect_true(all(is.finite(predict(f, f$design$data))))
  }
})

test_that("complete-workflow replication refits retention from baseline rosters", {
  skip_if_not_installed("survey")
  z <- attrition_fixture(); spec <- svy_attrition("observed", "y", analysis = "target")
  f <- svy_fit(z, "y", "x", family = "gaussian", v = 3, inner_v = 2, attrition = spec)
  ww <- cbind(z$w, z$w * ifelse(z$data$psu %% 2 == 0, .8, 1.2))
  rd <- survey::svrepdesign(weights = ~w, repweights = ww, data = z$data,
    type = "other", combined.weights = TRUE, scale = 1, rscales = c(.5, .5), mse = TRUE)
  a <- svy_refit_replicates(f, rd)
  expect_equal(a$status, rep("ok", 2))
  expect_equal(unname(a$estimates[1, ]), unname(f$metrics[1, ]))
  changed <- z; changed$w <- ww[, 2]; changed$data$w <- ww[, 2]
  oracle <- svy_fit(changed, "y", "x", family = "gaussian", v = 3, inner_v = 2, attrition = spec)
  expect_equal(unname(a$estimates[2, ]), unname(oracle$metrics[1, ]))
  expect_true(all(is.finite(a$variance)))
})

test_that("stacking and repeats retain their own cross-fitted retention weights", {
  z <- attrition_fixture(); z$data$binary <- as.numeric(z$data$y > 0)
  intercept <- svy_learner("custom", fit = function(x, y, w) weighted.mean(y, w),
    predict = function(model, x) rep(model, nrow(x)))
  f <- svy_fit(z, "binary", "x", learners = list(svy_learner(), intercept),
    v = 3, repeats = 2, inner_v = 2, stack = TRUE,
    attrition = svy_attrition("observed", "y", analysis = "target", stabilize = TRUE))
  expect_equal(rowSums(f$attrition$weights) / 2, f$design$w)
  expect_equal(rowMeans(f$predictions), f$consensus)
  expect_equal(sum(f$final$candidate_weights), 1)
  for (r in 1:2) expect_equal(f$metrics[r, ],
    svy_metrics(f$design$data$binary, f$predictions[, r], f$attrition$weights[, r]))
})

test_that("the retention helper accepts optional regularized learners", {
  skip_if_not_installed("glmnet")
  z <- attrition_fixture()
  a <- svy_attrition_weights(z, svy_attrition("observed", c("y", "x"), learner = svy_learner("glmnet")))
  expect_true(all(is.finite(a$weights)))
  expect_true(all(vapply(a$model$models, function(m) m$fit$spec$engine == "glmnet", logical(1))))
})
