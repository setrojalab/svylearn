test_that("tied weighted AUC agrees with pairwise probability", {
  y <- c(1,0,1,0,1,0); p <- c(.2,.2,.7,.7,.5,.1); w <- c(1,4,2,3,8,1)
  cmp <- outer(p[y == 1], p[y == 0], function(a,b) (a > b) + .5 * (a == b))
  truth <- sum(cmp * outer(w[y == 1], w[y == 0])) / (sum(w[y == 1]) * sum(w[y == 0]))
  expect_equal(unname(svy_metrics(y,p,w)["auc"]), truth, tolerance = 1e-14)
  expect_equal(svy_metrics(y,p,w), svy_metrics(y,p,w * 1e8), tolerance = 1e-14)
  expect_equal(unname(svy_metrics(y, rep(.5,6), w)["average_precision"]), sum(w*y)/sum(w))
  expect_true(is.na(svy_metrics(rep(0,6),p,w)["auc"]))
  expect_true(all(is.finite(svy_metrics(c(0,1),c(1,0))[1:2])))
  expect_error(svy_metrics(c(0,NA),c(.1,.5)), "nonmissing")
  expect_error(svy_metrics(y,p,-w), "Weights")
  expect_error(svy_metrics(y, p + 2,w), "Probabilities")
  expect_equal(svy_metrics(y,p,w * c(0,1,1,1,1,1)), svy_metrics(y[-1],p[-1],w[-1]))
})

test_that("whole clusters and panel subjects never cross folds; RNG is restored", {
  d <- svy_example(200, waves = 3); z <- svy_design(d,"weight","psu","stratum",subject="subject",wave="wave")
  set.seed(123); old <- .Random.seed; f <- svy_split(z,5,3,19)
  expect_identical(.Random.seed,old)
  expect_identical(f, svy_split(z,5,3,19))
  for (r in 1:3) {
    expect_true(all(vapply(split(f[,r],z$cluster), function(x) length(unique(x)) == 1, logical(1))))
    expect_true(all(vapply(split(f[,r],d$subject), function(x) length(unique(x)) == 1, logical(1))))
    expect_setequal(unique(f[,r]),1:5)
  }
  d$psu[2] <- 99
  expect_error(svy_design(d,"weight","psu","stratum",subject="subject"), "crosses")
  expect_error(svy_design(rbind(d,d[1,]),"weight",subject="subject",wave="wave"),"Duplicate")
})

test_that("preprocessing is frozen and novel categories do not change dimensions", {
  tr <- data.frame(x=c(1,2,NA,4), c=c("a","b",NA,"a"))
  rec <- svy_recipe(tr,c(1,1,5,1)); before <- serialize(rec,NULL)
  te <- data.frame(x=c(NA,100),c=c("new",NA))
  m <- bake_svy_recipe(rec,te)
  expect_identical(serialize(rec,NULL),before)
  expect_identical(colnames(m), colnames(bake_svy_recipe(rec,tr)))
  expect_true(all(is.finite(m)))
  expect_equal(rec$specs$x$median,2)
  expect_error(svy_recipe(data.frame(x=c(NA_real_,NA_real_))),"no observed")
})

test_that("outer evaluation outcomes cannot influence their fitted pipeline", {
  d <- svy_example(200); z <- svy_design(d,"weight","psu","stratum")
  f <- svy_fit(z,"y",c("x1","x2"),v=3,inner_v=2,seed=71)
  i <- which(f$folds[,1] == 1); d$y[i] <- 1 - d$y[i]
  f2 <- svy_fit(svy_design(d,"weight","psu","stratum"),"y",c("x1","x2"),v=3,inner_v=2,seed=71)
  expect_equal(f$predictions[i,1],f2$predictions[i,1],tolerance=1e-13)
  expect_equal(dim(f$predictions),c(200L,1L))
  expect_true(all(is.finite(predict(f,d))))
  expect_error(svy_fit(z,"y",c("x1","y")),"exclude")
})

test_that("calibration handles ties, monotonicity and population means", {
  y <- c(0,1,0,1,1,0); p <- c(.1,.2,.2,.3,.6,.9); w <- c(1,2,3,1,4,1)
  iso <- svy_calibrate(y,p,w,"isotonic",min_events=0)
  expect_true(all(diff(predict(iso,seq(0,1,length.out=100))) >= 0))
  expect_equal(predict(iso,c(.2,.2))[1],predict(iso,c(.2,.2))[2])
  a <- svy_align_prevalence(p,w,.17)
  expect_equal(sum(w*a)/sum(w),.17,tolerance=1e-10)
  const <- svy_calibrate(rep(0,6),p,w)
  expect_match(const$status,"fallback")
  expect_true(all(predict(const,p) > 0))
  expect_equal(predict(svy_calibrate(y,p,w,"none"),p),p)
})

test_that("zero and one-event folds return explicit finite fallback predictions", {
  d <- svy_example(200); d$y <- 0; d$y[1] <- 1
  z <- svy_design(d,"weight","psu","stratum")
  f <- svy_fit(z,"y",c("x1","x2"),learners=list(svy_learner("glmnet")),v=3,inner_v=2)
  expect_true(all(is.finite(f$predictions)))
  expect_true(all(f$predictions > 0 & f$predictions < 1))
  expect_true(all(grepl("fallback",f$audit$learner_status)))
  expect_true(all(grepl("fallback",f$audit$calibration_status)))
})

test_that("Taylor score estimates and SE match survey including domains", {
  skip_if_not_installed("survey")
  d <- svy_example(400); p <- plogis(-2 + d$x1); d$loss <- (d$y-p)^2
  z <- svy_design(d,"weight","psu","stratum")
  ds <- survey::svydesign(ids=~psu,strata=~stratum,weights=~weight,data=d,nest=TRUE)
  target <- survey::svymean(~loss,ds)
  got <- svy_infer(z,d$y,p)$scores
  expect_equal(got$estimate[2],unname(coef(target)),tolerance=1e-13)
  expect_equal(got$se[2],sqrt(unname(vcov(target)[1,1])),tolerance=1e-13)
  dom <- d$x1 > 0
  target <- survey::svymean(~loss,ds[dom,])
  got <- svy_infer(z,d$y,p,domain=dom)$scores
  expect_equal(got$se[2],sqrt(unname(vcov(target)[1,1])),tolerance=1e-13)
  rd <- svy_replicates(z,replicates=30)
  out <- svy_infer(z,d$y,p,replicate_design=rd)
  expect_equal(out$ranks$estimate,unname(svy_metrics(d$y,p,d$weight)[c("auc","average_precision")]))
  expect_true(all(out$ranks$se > 0))
})

test_that("optional engines, stacking and reporting run", {
  skip_if_not_installed("glmnet")
  d <- svy_example(600); z <- svy_design(d,"weight","psu","stratum")
  f <- svy_fit(z,"y",c("x1","x2"),learners=list(svy_learner(),svy_learner("glmnet")),
               v=3,inner_v=2,stack=TRUE)
  expect_equal(sum(f$final$candidate_weights),1,tolerance=1e-10)
  imp <- svy_importance(f,repeats=1)
  expect_equal(nrow(imp),2)
  expect_true(all(is.finite(imp$log_loss_increase)))
  for (fmt in c("csv","markdown","latex")) {
    path <- tempfile(); svy_export(f,path,fmt); expect_true(file.info(path)$size > 0)
  }
  for (engine in c("ranger","nnet")) {
    skip_if_not_installed(engine)
    spec <- svy_learner(engine,if(engine=="ranger") list(num.trees=20L) else list(size=2L,maxit=1000L))
    ff <- svy_fit(z,"y",c("x1","x2"),learners=list(spec),v=3,inner_v=2,calibration="none")
    expect_true(all(is.finite(predict(ff,d))))
  }
})

test_that("finite-population simulation reports MCSE and preserves seeds", {
  skip_if_not_installed("survey")
  d <- svy_example(200); d$p <- plogis(-2+d$x1)
  z <- svy_plasmode(d,"y","p","psu","stratum",replicates=4)
  expect_equal(nrow(z$replicates),8)
  expect_true(all(is.finite(z$summary$mcse_coverage)))
})
