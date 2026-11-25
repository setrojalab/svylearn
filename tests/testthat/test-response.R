response_data <- function() {
  d <- svy_example(240, seed = 29)
  set.seed(9)
  d$amount <- 2 + 1.5 * d$x1 - .6 * d$x2 + rnorm(nrow(d))
  d$count <- rpois(nrow(d), exp(.3 + .3 * d$x1))
  d$class <- factor(sample(c("z", "a", "q"), nrow(d), TRUE), levels = c("z", "a", "q"))
  d
}

test_that("response scores match independent algebra and Poisson family deviance", {
  y <- c(-2, 0, 4, 7); p <- c(-1, 1, 3, 5); w <- c(1, 3, 4, 2)
  got <- svy_metrics(y, p, w, "gaussian")
  expect_equal(unname(got["mse"]), sum(w * (y - p)^2) / sum(w))
  expect_equal(unname(got["mae"]), sum(w * abs(y - p)) / sum(w))
  expect_equal(unname(got["bias"]), sum(w * (p - y)) / sum(w))
  expect_equal(unname(got["rmse"]), sqrt(sum(w * (y - p)^2) / sum(w)))
  y <- c(0, 1, 4, 7); p <- c(.1, 1, 3, 5)
  expect_equal(unname(svy_metrics(y, p, w, "poisson")[1]),
    sum(stats::poisson()$dev.resids(y, p, w)) / sum(w))
  y <- factor(c("a", "b", "c", "a")); p <- rbind(c(.7,.2,.1),c(.1,.6,.3),c(.2,.2,.6),c(.4,.5,.1))
  colnames(p) <- levels(y)
  target <- diag(3)[as.integer(y), ]
  got <- svy_metrics(y, p[, c(3,1,2)], w, "multinomial")
  expect_equal(unname(got["log_loss"]), -sum(w * log(c(.7,.6,.6,.4))) / sum(w))
  expect_equal(unname(got["brier"]), sum(w * rowSums((target - p)^2)) / sum(w))
  expect_equal(unname(got["accuracy"]), .8)
  expect_equal(got, svy_metrics(y, p, w * 1e6, "multinomial"))
  expect_error(svy_metrics(c(1, -.5), c(1, 1), family = "poisson"), "integer")
  expect_error(svy_metrics(c(1, NA), c(1, 1), family = "gaussian"), "finite")
  expect_error(svy_metrics(as.character(y), p, family = "multinomial"), "factor")
  expect_error(svy_metrics(y, p * 2, family = "multinomial"), "row-normalized")
  expect_error(svy_metrics(y, unname(p), family = "multinomial"), "columns")
  expect_error(svy_metrics(1:2, c(1,-1), family = "poisson"), "predictions")
})

test_that("continuous and count final fits agree with direct weighted GLMs", {
  d <- response_data(); z <- svy_design(d, "weight", "psu", "stratum")
  for (family in c("gaussian", "poisson")) {
    response <- if (family == "gaussian") "amount" else "count"
    f <- svy_fit(z, response, c("x1", "x2"), family = family, calibration = "none", v = 3, inner_v = 2)
    ref <- stats::glm(stats::reformulate(c("x1", "x2"), response), data = d,
      weights = weight, family = if (family == "gaussian") stats::gaussian() else stats::poisson())
    expect_equal(predict(f, d), unname(predict(ref, d, type = "response")), tolerance = 1e-8)
    expect_equal(dim(f$predictions), c(nrow(d), 1L))
    expect_equal(summary(f)$family, family)
    expect_output(print(f), family)
    expect_true(all(is.finite(svy_importance(f, repeats = 1)[, 4])))
  }
})

test_that("multinomial GLM agrees with an independent nnet optimizer", {
  skip_if_not_installed("nnet")
  d <- response_data(); x <- as.matrix(d[c("x1","x2")])
  f <- .response_engine(x,d$class,d$weight,svy_learner("glm",list(lambda=0)),"multinomial")
  ref <- nnet::multinom(class~x1+x2,data=d,weights=d$weight/mean(d$weight),trace=FALSE,
    maxit=1000,reltol=1e-12)
  expect_equal(.response_predict_engine(f,x),unname(predict(ref,d,type="probs")) |>
    `colnames<-`(levels(d$class)),tolerance=1e-5)
  expect_error(.response_engine(x,d$class,d$weight,svy_learner("glm",list(lambda=-1)),"multinomial"),"penalty")
  expect_error(.response_engine(x,d$class,d$weight,svy_learner("glm",list(other=1)),"multinomial"),"accepts")
  # The default iteration budget must also handle nearly separated classes.
  iris_x <- bake_svy_recipe(svy_recipe(iris[1:4]),iris)
  near <- .response_engine(iris_x,iris$Species,rep(1,nrow(iris)),svy_learner(),"multinomial")
  expect_true(all(is.finite(.response_predict_engine(near,iris_x))))
})

test_that("family calibration obeys independent weighted moment identities", {
  y <- c(-1,2,5,3,7); p <- c(0,1,2,3,4); w <- c(1,2,4,1,2)
  cal <- .response_calibrator(y,p,w,"gaussian","auto")
  residual <- y-.response_calibrate_predict(cal,p)
  expect_equal(sum(w*residual),0,tolerance=1e-12)
  expect_equal(sum(w*residual*p),0,tolerance=1e-12)
  y <- c(0,1,3,2,6); p <- c(.2,.5,2,2,3)
  cal <- .response_calibrator(y,p,w,"poisson","auto")
  expect_equal(sum(w*.response_calibrate_predict(cal,p)),sum(w*y),tolerance=1e-12)
  y <- factor(c("a","b","c","a","c")); p <- rbind(c(.9,.05,.05),c(.8,.1,.1),
    c(.8,.1,.1),c(.5,.4,.1),c(.2,.3,.5)); colnames(p)<-levels(y)
  cal <- .response_calibrator(y,p,w,"multinomial","auto")
  pp <- .response_calibrate_predict(cal,p)
  expect_lte(svy_metrics(y,pp,w,"multinomial")[1],svy_metrics(y,p,w,"multinomial")[1]+1e-12)
  expect_identical(max.col(pp,ties.method="first"),max.col(p,ties.method="first"))
  expect_equal(rowSums(pp),rep(1,length(y)),tolerance=1e-12)
})

test_that("outer outcomes are isolated for every response family and stochastic learner", {
  d <- response_data()
  spec <- svy_learner("custom", fit = function(x,y,w) {
    if (is.factor(y)) list(p = setNames(rep(1 / nlevels(y), nlevels(y)), levels(y))) else
      list(p = sum(w * y) / sum(w) + runif(1))
  }, predict = function(model,x) if (length(model$p) > 1L)
    matrix(rep(model$p, each = nrow(x)), nrow(x), dimnames = list(NULL, names(model$p))) else rep(model$p,nrow(x)))
  for (family in c("gaussian", "poisson", "multinomial")) {
    response <- switch(family, gaussian = "amount", poisson = "count", multinomial = "class")
    z <- svy_design(d, "weight", "psu", "stratum")
    set.seed(83); before <- .Random.seed
    f <- svy_fit(z, response, c("x1", "x2"), learners = list(spec, svy_learner()),
      family = family, v = 3, inner_v = 2, repeats = 2, seed = 77, stack = TRUE)
    expect_identical(.Random.seed, before)
    i <- which(f$folds[, 1L] == 1L); dd <- d
    if (family == "multinomial") dd[[response]][i] <- levels(d[[response]])[(as.integer(d[[response]][i]) %% 3) + 1L] else
      dd[[response]][i] <- dd[[response]][i] + 10
    ff <- svy_fit(svy_design(dd, "weight", "psu", "stratum"), response, c("x1", "x2"),
      learners = list(spec, svy_learner()), family = family, v = 3, inner_v = 2, repeats = 2, seed = 77, stack = TRUE)
    if (family == "multinomial") {
      expect_equal(f$predictions[i, , 1], ff$predictions[i, , 1], tolerance = 1e-13)
      expect_equal(rowSums(f$consensus), rep(1, nrow(d)), tolerance = 1e-13)
      expect_equal(dim(f$predictions), c(nrow(d), 3L, 2L))
    } else expect_equal(f$predictions[i, 1], ff$predictions[i, 1], tolerance = 1e-13)
    expect_equal(predict(unserialize(serialize(f, NULL)), d), predict(f, d), tolerance = 0)
    expect_equal(sum(f$final$candidate_weights), 1, tolerance = 1e-13)
  }
})

test_that("all documented engine-family combinations predict valid responses", {
  d <- response_data(); x <- as.matrix(d[c("x1", "x2")])
  for (engine in c("glm", "glmnet", "ranger", "xgboost", "nnet")) {
    if (engine != "glm" && !requireNamespace(engine, quietly = TRUE)) next
    parameters <- switch(engine, ranger = list(num.trees = 20), xgboost = list(nrounds = 5),
      nnet = list(size = 1, decay = .1, maxit = 2000), list())
    for (family in c("gaussian", "poisson", "multinomial")) {
      if (engine == "nnet" && family == "poisson") next
      y <- d[[switch(family, gaussian = "amount", poisson = "count", multinomial = "class")]]
      f <- .response_engine(x, y, d$weight, svy_learner(engine, parameters), family)
      p <- .response_predict_engine(f, x)
      expect_true(all(is.finite(p)), info = paste(engine, family))
      expect_equal(p, .response_predict_engine(unserialize(serialize(f, NULL)), x), tolerance = 0)
      if (family == "multinomial") {
        expect_equal(colnames(p), levels(y))
        expect_equal(rowSums(p), rep(1, length(y)), tolerance = 1e-6)
        expect_equal(dim(.response_predict_engine(f, x[1, , drop = FALSE])), c(1L, 3L))
      } else expect_length(p, length(y))
    }
  }
})

test_that("sparse response fallbacks and failed candidates are explicit", {
  d <- response_data(); z <- svy_design(d, "weight", "psu", "stratum")
  bad <- svy_learner("custom", fit = function(x,y,w) stop("intentional failure"), predict = function(model,x) 0)
  f <- svy_fit(z, "amount", c("x1", "x2"), learners = list(bad, svy_learner()), family = "gaussian", v = 3, inner_v = 2)
  expect_true(all(is.infinite(f$audit$inner_loss[f$audit$candidate == 1])))
  expect_equal(f$final$candidate_weights, c(0,1))
  expect_error(svy_fit(z,"amount",c("x1","x2"),learners=list(bad),family="gaussian",v=3,inner_v=2), "All learner")
  for (family in c("gaussian", "poisson", "multinomial")) {
    dd <- d
    dd$constant <- if (family == "multinomial") factor(rep("a",nrow(d)), levels=c("a","b","c")) else rep(0,nrow(d))
    f <- svy_fit(svy_design(dd,"weight","psu","stratum"), "constant", c("x1","x2"), family=family,v=3,inner_v=2)
    expect_true(all(is.finite(predict(f,dd))))
    expect_true(all(grepl("fallback",f$audit$learner_status)))
    if (family == "multinomial") expect_true(all(predict(f,dd) > 0)) else expect_equal(predict(f,dd),rep(0,nrow(d)))
  }
  expect_error(svy_fit(z,"amount",c("x1","weight"),family="gaussian"),"Predictors")
  expect_error(svy_fit(z,"amount","x1",family="gaussian",inner_v=1),"inner_v")
  expect_error(svy_fit(z,"amount","x1",family="gaussian",learners=list()),"learners")
})

test_that("response-loss inference matches parent-domain and replicate survey references", {
  skip_if_not_installed("survey")
  d <- response_data(); z <- svy_design(d,"weight","psu","stratum")
  p <- 2 + d$x1; d$loss <- (d$amount - p)^2; domain <- d$x2 > 0
  ds <- survey::svydesign(ids=~psu,strata=~stratum,weights=~weight,nest=TRUE,data=d)
  ref <- survey::svymean(~loss, ds[domain,])
  got <- svy_infer(z,d$amount,p,domain=domain,family="gaussian")$scores
  expect_equal(got$estimate[1],unname(coef(ref)[1]),tolerance=1e-12)
  expect_equal(got$se[1],sqrt(unname(vcov(ref)[1,1])),tolerance=1e-12)
  rd <- svy_replicates(z,type="JKn"); rd$variables$loss <- d$loss
  ref <- survey::svymean(~loss,rd)
  got <- svy_infer(z,d$amount,p,replicate_design=rd,family="gaussian")$scores
  expect_equal(got$se[1],sqrt(unname(vcov(ref)[1,1])),tolerance=1e-12)
  for (family in c("poisson","multinomial")) {
    y <- if(family=="poisson") d$count else d$class
    p <- if(family=="poisson") exp(.3+.3*d$x1) else matrix(1/3,nrow(d),3,dimnames=list(NULL,levels(y)))
    got <- svy_infer(z,y,p,family=family)
    expect_equal(got$scores$estimate,unname(svy_metrics(y,p,d$weight,family)[seq_len(nrow(got$scores))]))
  }
  f <- svy_fit(z,"amount",c("x1","x2"),family="gaussian",v=3,inner_v=2)
  rr <- svy_refit_replicates(f,svy_replicates(z,type="bootstrap",replicates=5,seed=81))
  expect_true(all(rr$status=="ok")); expect_true(all(is.finite(rr$variance)))
  expect_error(svy_infer(z,d$amount,2+d$x1,domain=rep(FALSE,nrow(d)),family="gaussian"),"domain")
  expect_error(svy_infer(z,d$amount,2+d$x1,level=2,family="gaussian"),"confidence")
  audit <- svy_subgroups(z,d$amount,2+d$x1,domain,family="gaussian")
  expect_true(all(vapply(audit,function(x)x$sufficient,logical(1))))
  expect_equal(audit$`TRUE`$inference$scores$se,
    svy_infer(z,d$amount,2+d$x1,domain=domain,family="gaussian")$scores$se)
  probs <- matrix(1/3,nrow(d),3,dimnames=list(NULL,levels(d$class)))
  audit <- svy_subgroups(z,d$class,probs,d$class,family="multinomial")
  expect_true(all(vapply(audit,function(x)!x$sufficient && is.null(x$inference),logical(1))))
  expect_error(svy_subgroups(z,d$amount,2+d$x1,domain,min_events=0,family="gaussian"),"threshold")
})
