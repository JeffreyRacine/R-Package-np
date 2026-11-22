test_that("regression constructions preserve finite-grid scope and caller state", {
  old <- options(np.messages=FALSE, np.plot.progress=FALSE)
  on.exit(options(old), add=TRUE)
  x <- data.frame(x=seq(-1,1,length.out=48))
  y <- x$x^2 + .15*sin(seq_len(nrow(x)))
  b <- npregbw(xdat=x,ydat=y,bws=.35,regtype="ll",bandwidth.compute=FALSE)
  fit <- npreg(b,txdat=x,tydat=y)
  set.seed(20538);seed <- .Random.seed
  run <- function(...) suppressWarnings(plot(fit,errors="bootstrap",neval=4L,
    B=59L,random.seed=231L,output="data",...))
  omit <- run(); normal <- run(interval="normal-sd")
  expect_identical(omit$panels,normal$panels)
  expect_identical(normal$request$band,"pointwise")
  expect_false("interval" %in% omit$request$supplied)
  expect_true("interval" %in% normal$request$supplied)
  null <- run(interval=NULL)
  expect_true("interval" %in% null$request$supplied)
  expect_identical(null$panels,omit$panels)
  joint <- run(interval="percentile",band="all")
  expect_setequal(names(joint$inference[[1L]]$components),
    c("pointwise","simultaneous","bonferroni"))
  expect_identical(joint$inference[[1L]]$components$simultaneous$interval,"percentile")
  expect_error(run(band="pmzsd"),"interval='normal-sd'",fixed=TRUE)
  expect_error(run(interval="normal-sd",band="simultaneous"),"unavailable")
  expect_error(run(interval="asymmetric-sd"),"unavailable")
  expect_error(run(interval="symmetric-sd",boot.control=historical_boot_control(calibration="absolute")),"conflicting")
  expect_error(run(boot.control=historical_boot_control(calibration="default")),"omit interval")
  expect_error(run(plot.errors.boot.calibration="absolute"),"replace plot.errors.boot.calibration with interval",fixed=TRUE)
  expect_error(run(interval=c("percentile","normal-sd")),"single supported")
  expect_error(run(interval="percentile",bias.control=1),"np_bias_control",fixed=TRUE)
  expect_identical(.Random.seed,seed)
})

test_that("regression residual percentile components identify their actual ranks", {
  old <- options(np.messages=FALSE, np.plot.progress=FALSE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(x=seq(-1,1,length.out=44L))
  y <- sin(2*x$x)+.25*cos(seq_len(nrow(x)))
  b <- npregbw(xdat=x,ydat=y,bws=.45,regtype="ll",bandwidth.compute=FALSE)
  fit <- npreg(b,txdat=x,tydat=y)
  for (grid in c(1L,5L)) for (slope in c(FALSE,TRUE)) {
    z <- suppressWarnings(plot(fit,center="bias-corrected",errors="bootstrap",
      interval="outer-calibrated-percentile",band="all",gradients=slope,
      neval=grid,B=199L,random.seed=20552L,output="data"))
    ci <- z$inference[[1L]]$components
    expect_true(all(vapply(ci,function(v)identical(v$interval,"outer-calibrated-percentile"),logical(1))))
    if (grid==1L) {
      expect_match(ci$pointwise$endpoint_construction,"single-target")
      for(key in c("pointwise","bonferroni")){
        a <- as.data.frame(z,band=key); s <- as.data.frame(z,band="simultaneous")
        expect_identical(a$lower,s$lower);expect_identical(a$upper,s$upper)
      }
    } else {
      expect_identical(ci$pointwise$calibration,if(slope)"E6 outer rank calibration" else "nominal marginal ranks")
      expect_match(ci$simultaneous$endpoint_construction,if(slope)"E6" else "E1")
      expect_match(ci$bonferroni$endpoint_construction,"Bonferroni")
    }
  }
})
