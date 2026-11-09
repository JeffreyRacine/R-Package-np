test_that("stored fit owners refuse replacement subset without evaluating it", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
 n<-24L
  d<-data.frame(x=seq(.02,.98,length.out=n),z=cos(seq_len(n)*.73))
  d$y<-sin(4*d$x)+d$z/7
  d$fy<-factor(rep(c('a','b'),12))
  h<-.3
  makers <- list(
    npreg=function()npregbw(y~x,data=d,bws=h,bandwidth.compute=FALSE),
    npudens=function()npudensbw(~x,data=d,bws=h,bandwidth.compute=FALSE),
    npudist=function()npudistbw(~x,data=d,bws=h,bandwidth.compute=FALSE),
    npcdens=function()npcdensbw(y~x,data=d,bws=c(h,h),bandwidth.compute=FALSE),
    npcdist=function()npcdistbw(y~x,data=d,bws=c(h,h),bandwidth.compute=FALSE),
    npindex=function()npindexbw(y~x+z,data=d,bws=c(1,.5,h),bandwidth.compute=FALSE),
    npplreg=function()npplregbw(y~z|x,data=d,bws=matrix(c(h,h),2,1),bandwidth.compute=FALSE),
    npscoef=function()npscoefbw(y~z|x,data=d,bws=h,bandwidth.compute=FALSE),
    npconmode=function()npcdensbw(fy~x,data=d,bws=c(.2,h),bandwidth.compute=FALSE),
    npsigtest=function()npregbw(y~x,data=d,bws=h,bandwidth.compute=FALSE),
    npcopula=function()npudistbw(~x+z,data=d,bws=c(h,h),bandwidth.compute=FALSE),
    npreghat=function()npregbw(y~x,data=d,bws=h,bandwidth.compute=FALSE),
    npqreg=function()npcdistbw(y~x,data=d,bws=c(h,h),bandwidth.compute=FALSE))

  calls <- new.env(parent = emptyenv())
  for (family in names(makers)) {
    bws <- makers[[family]]()
    args <- list(bws = bws, data = d)
    if (family == "npcopula") {
      args$data <- d[c("x", "z")]
      args$u <- data.frame(x = c(.3, .6), z = c(.4, .7))
    }
    if (family == "npsigtest") args$B <- 9L
    if (family == "npreghat") args$output <- "apply"
    for (value in list(NULL, TRUE, d$x > .4,
                       quote({calls$n <- calls$n + 1L; stop("forced subset")}))) {
      calls$n <- 0L
      args["subset"] <- list(value)
      invocation <- as.call(c(list(as.name(family)), args))
      expect_error(eval(invocation), "replacement subset", fixed = TRUE,
                   info = family)
      expect_identical(calls$n, 0L, info = family)
    }
  }
  # Native conditional-mode dispatch must reach the same admission check.
  bws <- makers$npconmode()
  expect_error(npconmode(bws, txdat = d["x"], tydat = d["fy"],
                         subset = stop("forced subset")), "replacement subset")
})

test_that("constructor subsets, retained samples and S3 extensions keep ownership", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  d <- data.frame(x = seq(.02, .98, length.out = 24L))
  d$y <- sin(4 * d$x)
  calls <- new.env(parent = emptyenv()); calls$n <- 0L
  fit <- npreg(y ~ x, data = d, subset = {calls$n <- calls$n + 1L; x > .4},
               bws = .3, bandwidth.compute = FALSE)
  keep <- d$x > .4
  H <- dnorm(outer(d$x[keep], d$x[keep], "-") / .3)
  H <- H / rowSums(H)
  expect_identical(calls$n, 1L)
  expect_equal(as.numeric(fitted(fit)), drop(H %*% d$y[keep]), tolerance = 1e-12)
  expect_identical(fitted(npreg(fit$bws)), fitted(fit))
  native <- npregbw(xdat = d["x"], ydat = d$y, bws = .3,
                     bandwidth.compute = FALSE)
  fitted.native <- npreg(native)
  expect_error(npreg(native, subset = stop("forced")), "replacement subset")
  expect_error(npreghat(fitted.native, subset = stop("forced")), "replacement subset")
  expect_error(npsigtest(fitted.native, subset = stop("forced")), "replacement subset")
  custom <- native
  custom[["call"]] <- custom[[".np.native.training"]] <- NULL
  class(custom) <- "subset_test_owner"
  npreg.subset_test_owner <- function(bws, ..., subset) subset
  expect_identical(npreg(custom, txdat = d["x"], subset = "custom policy"),
                   "custom policy")
  lsq <- nplsqregbw(bws = y ~ x, data = d, bw = .3, bandwidth.compute = FALSE)
  expect_identical(fitted(nplsqreg(lsq, subset = NULL)), fitted(nplsqreg(lsq)))
  expect_identical(fitted(nplsqreg(lsq, subset = TRUE)), fitted(nplsqreg(lsq)))
  expect_error(nplsqreg(lsq, subset = d$x > .4), "retains its training sample")
  reg <- npregbw(y ~ x, data = d, bws = .3, bandwidth.compute = FALSE)
  expect_length(fitted(nplsqreg(reg, data = d, subset = x > .4)), sum(keep))
})
