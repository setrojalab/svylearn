test_that("weighted PCA reconstructs covariance and stays frozen for high-dimensional inputs", {
  set.seed(913); x <- as.data.frame(matrix(rnorm(120),40,3)); w <- seq_len(40)
  rec <- svy_recipe(x,w,scale=FALSE,pca_rank=3)
  centered <- sweep(as.matrix(x),2,colSums(as.matrix(x)*w)/sum(w))
  covariance <- crossprod(centered*sqrt(w/sum(w)))
  rebuilt <- rec$pca$rotation %*% diag(rec$pca$variance) %*% t(rec$pca$rotation)
  expect_equal(unname(rebuilt),unname(covariance),tolerance=1e-10)
  z <- bake_svy_recipe(rec,x)
  expect_equal(unname(colSums(z*w)/sum(w)),rep(0,3),tolerance=1e-12)
  expect_equal(unname(crossprod(z*sqrt(w/sum(w)))),diag(rec$pca$variance),tolerance=1e-10)
  wide <- as.data.frame(matrix(rnorm(20*200),20,200)); fit <- svy_recipe(wide,pca_rank=7)
  expect_equal(dim(bake_svy_recipe(fit,wide)),c(20L,7L))
  before <- fit$pca
  bake_svy_recipe(fit,wide*1000)
  expect_identical(fit$pca,before)
  expect_error(svy_recipe(x,pca_rank=0),"positive")
  expect_error(svy_recipe(data.frame(x=rep(1,5)),pca_rank=1),"nonconstant")
})

test_that("dependency groups merge transitively while preserving inference metadata", {
  d <- data.frame(w=1,id=1:12,h=rep(1:3,each=4),site=c(1,1,2,3,1,4,5,5,6,7,8,8),
    graph=c(1,2,3,3,4,5,6,7,8,9,10,11))
  a <- svy_design(d,"w","id","h",fold_groups=c("site","graph"))
  split <- svy_split(a,v=3,repeats=5)
  for(r in 1:5)for(v in c("site","graph"))
    expect_true(all(vapply(split(split[,r],d[[v]]),function(z)length(unique(z))==1L,logical(1))))
  expect_identical(a$ids,"id"); expect_identical(a$strata,"h")
  expect_equal(length(unique(a$cluster)),12L)
  expect_lt(length(unique(a$fold_cluster)),12L)
  b <- svylearn:::.subset_design(a,1:8)
  expect_equal(length(b$fold_cluster),8L)
  impossible <- svy_design(transform(d,site=1),"w","id","h",fold_groups="site")
  expect_error(svy_split(impossible,2),"Fewer independent")
})

test_that("trajectory features reject future information and duplicate waves", {
  d <- data.frame(id=rep(c("a","b","c"),each=3),time=rep(1:3,3),
    x=c(2,4,1000,1,NA,900,NA,NA,500))
  a <- svy_trajectory(d,"id","time","x",cutoff=2)
  expect_equal(a$x_n,c(2,1,0)); expect_equal(a$x_mean,c(3,1,NA))
  expect_equal(a$x_last,c(4,1,NA)); expect_equal(a$x_slope,c(2,NA,NA))
  d$x[d$time>2] <- -1e9
  expect_identical(svy_trajectory(d,"id","time","x",cutoff=2),a)
  expect_error(svy_trajectory(rbind(d,d[1,]),"id","time","x",2),"unique")
  expect_error(svy_trajectory(d,"id","time","x",Inf),"finite prespecified")
})

test_that("agency-specified replicate degrees of freedom are honored", {
  skip_if_not_installed("survey")
  d <- svy_example(200); z <- svy_design(d,"weight","psu","stratum")
  rd <- svy_replicates(z,type="bootstrap",replicates=100)
  wz <- svy_design(d,"weight","psu","stratum",repweights=as.matrix(weights(rd,type="analysis")),
    scale=rd$scale,rscales=rd$rscales,mse=rd$mse,replicate_df=80)
  expect_equal(survey::degf(as_survey_design(wz)),80)
  expect_error(svy_design(d,"weight",replicate_df=80),"requires supplied")
  expect_equal(svy_infer(wz,d$y,plogis(d$x1))$scores$df,rep(80,2))
})
