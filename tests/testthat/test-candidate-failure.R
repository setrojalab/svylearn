test_that("failed candidates are excluded and their failure is recorded", {
  d<-svy_example(200);z<-svy_design(d,"weight","psu","stratum")
  bad<-svy_learner("custom",fit=function(x,y,w)stop("declared test failure"),
                  predict=function(model,x)rep(.5,nrow(x)))
  for(stack in c(FALSE,TRUE)) {
    f<-svy_fit(z,"y",c("x1","x2"),learners=list(bad,svy_learner()),
      v=3,inner_v=2,stack=stack)
    expect_true(all(f$audit$mixture_weight[f$audit$candidate==1]==0))
    expect_true(all(grepl("declared test failure",f$audit$inner_status[f$audit$candidate==1])))
    expect_true(all(is.finite(f$predictions)))
  }
  expect_error(svy_fit(z,"y",c("x1","x2"),learners=list(bad),v=3,inner_v=2),"All learner candidates failed")
})
test_that("outer isolation includes outcome-dependent random-number consumption", {
  d<-svy_example(400)
  design<-svy_design(d,"weight","psu","stratum")
  stochastic<-svy_learner("custom",
    fit=function(x,y,w) {runif(1+sum(y));runif(ncol(x),-.5,.5)},
    predict=function(model,x)plogis(as.numeric(x%*%model)-1))
  a<-svy_fit(design,"y",c("x1","x2"),list(stochastic),v=4,inner_v=3,seed=18)
  held<-which(a$folds[,1]==3L)
  altered<-design;altered$data$y[held]<-1L-altered$data$y[held]
  b<-svy_fit(altered,"y",c("x1","x2"),list(stochastic),v=4,inner_v=3,seed=18)
  expect_equal(a$predictions[held,1],b$predictions[held,1],tolerance=1e-14)
})
