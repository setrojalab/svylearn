test_that("multistage FPC, supplied replicates and raking match survey", {
  skip_if_not_installed("survey")
  d <- svy_example(400); d$fpc1 <- 20; d$fpc2 <- 50
  d$unit <- ave(d$subject,interaction(d$stratum,d$psu),FUN=seq_along)
  d$p <- plogis(-2+d$x1); d$loss <- (d$y-d$p)^2
  z <- svy_design(d,"weight",c("psu","unit"),"stratum",c("fpc1","fpc2"))
  ds <- survey::svydesign(ids=~psu+unit,strata=~stratum,weights=~weight,fpc=~fpc1+fpc2,nest=TRUE,data=d)
  got <- svy_infer(z,d$y,d$p)$scores
  target <- survey::svymean(~loss,ds)
  expect_equal(got$se[2],sqrt(unname(vcov(target)[1,1])),tolerance=1e-12)
  z <- svy_design(d,"weight","psu","stratum")
  rd <- svy_replicates(z,type="JKn")
  zr <- svy_design(d,"weight","psu","stratum",repweights=as.matrix(weights(rd,type="analysis")),
    scale=rd$scale,rscales=rd$rscales,mse=rd$mse)
  target <- survey::svymean(~loss,rd)
  expect_equal(svy_infer(zr,d$y,d$p)$scores$se[2],sqrt(unname(vcov(target)[1,1])),tolerance=1e-12)
  calibrated <- svy_rake(z,~x1,population=c(`(Intercept)`=sum(d$weight),x1=0))
  expect_equal(unname(coef(survey::svytotal(~x1,calibrated))),0,tolerance=1e-10)
})

test_that("Fay, BRR and bootstrap replication preserve fixed points", {
  skip_if_not_installed("survey")
  d <- svy_example(400); d <- d[d$psu<=2,]
  z <- svy_design(d,"weight","psu","stratum")
  p <- plogis(-2+d$x1)
  for (scheme in c("BRR","Fay","bootstrap")) {
    rd <- if(scheme=="Fay") svy_replicates(z,type=scheme,fay.rho=.3) else svy_replicates(z,type=scheme,replicates=30)
    expect_equal(svy_infer(z,d$y,p,replicate_design=rd)$ranks$estimate,
      unname(svy_metrics(d$y,p,d$weight)[c("auc","average_precision")]),tolerance=1e-12)
  }
})

test_that("singleton policy is explicit and options are restored", {
  skip_if_not_installed("survey")
  d <- svy_example(200); d <- d[!(d$stratum==1 & d$psu!=1),]
  z <- svy_design(d,"weight","psu","stratum")
  expect_equal(svy_diagnostics(z)$singleton_strata,1)
  old <- getOption("survey.lonely.psu")
  expect_error(svy_infer(z,d$y,plogis(d$x1)),"only one PSU")
  expect_identical(getOption("survey.lonely.psu"),old)
  z$lonely.psu <- "adjust"
  expect_true(all(is.finite(svy_infer(z,d$y,plogis(d$x1))$scores$se)))
  expect_identical(getOption("survey.lonely.psu"),old)
})

test_that("curves integrate to the tie-corrected AUC", {
  y<-c(0,1,1,0,1,0);p<-c(.1,.2,.2,.6,.9,.9);w<-c(1,2,3,4,2,1)
  expect_equal(svy_curves(y,p,w,max_fpr=1)$partial_auc,
               unname(svy_metrics(y,p,w)["auc"]),tolerance=1e-14)
  expect_equal(svy_curves(y,rep(.5,6),w,max_fpr=.1)$partial_auc,.005,tolerance=1e-14)
  expect_error(svy_curves(rep(0,6),p,w),"Both")
  ce <- svy_calibration_error(y,p,w)
  expect_true(is.finite(ce$ICI))
})

test_that("complete workflow refits respect replicate conventions", {
  skip_if_not_installed("survey")
  d<-svy_example(300);z<-svy_design(d,"weight","psu","stratum")
  f<-svy_fit(z,"y",c("x1","x2"),v=3,inner_v=2)
  rd<-svy_replicates(z,replicates=4)
  rr<-svy_refit_replicates(f,rd)
  expect_true(all(rr$status=="ok"))
  expect_equal(dim(rr$variance),c(4L,4L))
})

test_that("XGBoost adapter works in strict mode", {
  skip_if_not_installed("xgboost")
  old<-options(xgboost.strict_mode=TRUE);on.exit(options(old))
  d<-svy_example(300)
  f<-svy_fit(svy_design(d,"weight","psu","stratum"),"y",c("x1","x2"),
    learners=list(svy_learner("xgboost",list(nrounds=5))),v=3,inner_v=2)
  expect_true(all(is.finite(predict(f,d))))
  expect_type(f$final$fits[[1]]$fit$model,"raw")
  path<-tempfile();saveRDS(f,path)
  expect_equal(predict(readRDS(path),d),predict(f,d),tolerance=1e-12)
})
