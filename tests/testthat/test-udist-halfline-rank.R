test_that("boundary-rank uses the known support and preserves the signed estimate", {
  x <- data.frame(x = seq(.08, 4, length.out=43L)^1.3)
  bw <- npudistbw(dat=x, bws=.38, bandwidth.compute=FALSE,
                  ckerbound="fixed", ckerlb=0, ckerub=Inf)
  bc <- np_udist_bias_control(calibration="boundary-rank")
  expect_error(np_udist_bias_control("cdf", "boundary-rank"), "joint-process")
  prep <- function(b=bw, d=x, control=bc, errors="asymptotic", band="all")
    .np_udist_plot_prepare(control,b,d,"bias-corrected",errors,"inid",band)
  a <- prep()
  expect_identical(a$metadata$kernel, "gaussian2-halfline")
  expect_error(prep(control=NULL),"unbounded")
  expect_error(prep(errors="bootstrap"),"boundary-rank")
  expect_error(prep(band="pmzsd"),"boundary-rank")
  expect_error(prep(d=x[c(1L,seq_len(nrow(x)-1L)),,drop=FALSE]),"untied")
  g <- data.frame(x=c(0, .001, min(x$x), .35, 1, max(x$x)))
  r <- .np_udist_bc_rows(a,g)
  b <- npudistbw(dat=data.frame(x=7-2*x$x), bws=.76, bandwidth.compute=FALSE,
                 ckerbound="fixed", ckerlb=-Inf, ckerub=7)
  up <- .np_udist_bc_rows(prep(b,data.frame(x=7-2*x$x)),data.frame(x=7-2*g$x))
  native.upper <- fitted(npudist(bws=b, tdat=data.frame(x=7-2*x$x),
                                edat=data.frame(x=7-2*g$x)))
  expect_equal(up$original.center, as.numeric(native.upper), tolerance=2e-8)
  expect_equal(r$se[4L], sd(r$corrected[4L,])/sqrt(nrow(x)), tolerance=1e-13)
  expect_equal(r$center, 1-up$center, tolerance=1e-12)
  expect_equal(r$se, up$se, tolerance=1e-12)
  expect_equal(r$factor,-up$factor, tolerance=1e-12)
  expect_equal(r$original.center, fitted(npudist(bws=bw,tdat=x,edat=g)), tolerance=2e-8, ignore_attr=TRUE)
  expect_identical(as.numeric(r$center[1]),0)
  expect_identical(as.numeric(r$se[1]),0)
  expect_error(.np_udist_bc_rows(a,data.frame(x=-.01)), "within")
  inf <- .np_udist_joint_inference(a,g,B=99L,alpha=c(.05,.1))
  for (i in 1:2) {
    v <- inf$intervals
    for (s in names(v)) {
      expect_equal(v[[s]][1,,i],c(lower=0,upper=0))
      expect_true(all(v[[s]][,1,i]<=r$center & r$center<=v[[s]][,2,i]))
    }
    expect_true(all(v$bonferroni[,1,i]<=v$simultaneous[,1,i] & v$simultaneous[,1,i]<=v$pointwise[,1,i]))
    expect_true(all(v$bonferroni[,2,i]>=v$simultaneous[,2,i] & v$simultaneous[,2,i]>=v$pointwise[,2,i]))
  }
})

test_that("boundary-rank includes the empirical bound for infeasible signed intervals", {
  x <- c(.1, .2, .5, .7, .9)
  grid <- c(0, .1, .3, .6, .9)
  center <- c(0, -.02, -.4, 1.4, 1.03)
  alpha <- c(.05, .10)
  point <- matrix(.01, length(grid), length(alpha))
  base <- array(NA_real_, c(length(grid), 2L, length(alpha)))
  for (j in seq_along(alpha)) {
    base[, 1L, j] <- center-point[, j]
    base[, 2L, j] <- center+point[, j]
  }
  answer <- list(center=center, intervals=list(pointwise=base),
                 critical=list(pointwise=point), settings=list(B=99L))
  got <- .np_udist_halfline_envelope(answer, list(pointwise=point),
                                     x, grid, alpha, 1L)
  for (j in seq_along(alpha)) {
    v <- got$intervals$pointwise[, , j]
    radius <- sqrt(log(2/alpha[j])/(2*length(x)))
    F <- vapply(grid[3:4], function(g) mean(x <= g), 0.0)
    expect_true(all(v[3:4,1L] <= pmax(0, F-radius)))
    expect_true(all(v[3:4,2L] >= pmin(1, F+radius)))
    expect_equal(v[1L,], c(0,0))
    expect_gte(v[2L,2L], stats::qbeta(1-alpha[j]/2,1,length(x)))
    expect_lte(v[5L,1L], stats::qbeta(alpha[j]/2,length(x),1))
    expect_lte(v[2L,1L], stats::qbeta(alpha[j]/2,1,length(x)))
    expect_gte(v[5L,2L], stats::qbeta(1-alpha[j]/2,length(x),1))
    expect_equal(got$critical$pointwise[, j], pmax(center-v[,1L], v[,2L]-center))
  }
  expect_identical(got$center, center)
})

test_that("half-line queries outside the sample use the common extrema events", {
  x <- c(.1,.2,.5,.7,.9)
  g <- c(.03,.08,.95)
  center <- c(-.002,.004,1.01)
  width <- matrix(.001,3L,1L)
  v <- array(c(center-.001, center+.001), c(3L,2L,1L))
  raw <- list(center=center, intervals=list(pointwise=v), critical=list(pointwise=width), settings=list())
  ans <- .np_udist_halfline_envelope(raw,list(pointwise=width),x,g,.05,integer())
  got <- ans$intervals$pointwise[, , 1L]
  expect_equal(got[1:2,2L], rep(stats::qbeta(.975,1L,5L),2))
  expect_equal(got[3L,1L],stats::qbeta(.025,5L,1L))
  expect_identical(ans$settings$critical.components[[1L]]$rank.upper,1:2)
  expect_identical(ans$settings$critical.components[[1L]]$feasibility,integer())
})
