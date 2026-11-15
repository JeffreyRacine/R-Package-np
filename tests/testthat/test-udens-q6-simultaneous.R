test_that("Q6 compact operators retain their reproducing moments", {
  for (kernel in c("uniform", "epanechnikov")) {
    span <- if (kernel == "uniform") 1 else sqrt(5)
    for (k in 0:5) {
      value <- integrate(function(r) r^k * .np_udens_q6_kernel(r, kernel),
                         -span, span, abs.tol = 1e-12, rel.tol = 1e-11)$value
      expect_equal(value, as.numeric(k == 0), tolerance = 1e-11)
    }
  }
})

test_that("Q6 uses the complete Student6 reference and original scale", {
  x <- qnorm((1:20)/21)
  h <- .6
  grid <- c(-.7, .2, .8)
  for (kernel in c("uniform", "epanechnikov")) {
    covariance <- .np_udens_q6_covariance(x, h, grid, kernel, 96L)
    v <- mean((x-mean(x))^2)
    expect_equal(covariance$m4, mean((x-mean(x))^4) + 6*h*h*v + 6*h^4)
    ref <- .np_udens_q6_reference(x, h, grid, kernel)
    expect_equal(ref$se, sqrt(diag(covariance$C)), tolerance = 1e-12)
    expect_true(all(ref$diagnostics <= 1))
    # Independently integrate the reference density against the compact operator.
    span <- if (kernel == "uniform") 1 else sqrt(5)
    q <- function(u) vapply(u, function(t) mean(dt((t-x)/(h*sqrt(2/3)),6))/(h*sqrt(2/3)), 0.0)
    oracle <- vapply(grid, function(t) {
      integrate(function(r) .np_udens_q6_kernel(r,kernel) *
        (q(t+h*r)-dnorm(t+h*r,mean(x),sqrt(v+h*h))), -span, span,
        rel.tol=1e-11, abs.tol=1e-12)$value + dnorm(t,mean(x),sqrt(v+h*h))
    }, 0.0)
    expect_equal(ref$reference.center, oracle, tolerance = 1e-10)
    expect_false(isTRUE(all.equal(ref$center, ref$reference.center, tolerance=1e-9)))
    values <- c(.5, 1, 2, 3)
    out <- .np_udens_q6_intervals(ref, values, .1)
    critical <- unname(quantile(values,.9,type=7))
    expect_identical(out$critical, critical)
    expect_identical(out$lower, ref$center-critical*ref$se)
    expect_identical(out$upper, ref$center+critical*ref$se)
    perm <- .np_udens_q6_reference(x,h,grid[c(3,1,3,2)],kernel)
    expect_identical(perm$covariance, ref$covariance)
    expect_identical(.np_udens_q6_intervals(perm,values,.1)$lower,
                     out$lower[c(3,1,3,2)])
  }
})

test_that("Q6 numerical failures remain failures", {
  expect_error(.np_udens_q6_geometry(1:9,.1,0,"uniform"), "invalid Q6")
  expect_error(.np_udens_q6_geometry(rep(1,10),.1,0,"uniform"), "guide variance")
  expect_error(.np_udens_q6_geometry(1:10,0,0,"uniform"), "invalid Q6")
  expect_error(.np_udens_q6_budget(1000,100000,9999), "workspace bound")
  expect_error(.np_udens_q6_check(1,1.001,"test"), "refinement failed")
  expect_error(.np_udens_q6_check(Inf,1,"test"), "nonfinite")
  expect_error(.np_udens_q6_kernel(1,"gaussian"), "uniform or Epanechnikov")
  expect_error(.np_udens_q6_intervals(list(),c(1,Inf),.05), "maximum calibration")
})


test_that("Q6 blocks retain donor-first then noise stream and a global maximum", {
  .np_with_seed(42L, {
    x <- qnorm((seq_len(20L)-.5)/20); n <- length(x); B <- 17L; h <- .6
    start <- get(".Random.seed", .GlobalEnv)
    donor <- sample.int(n, n*B, TRUE); noise <- rnorm(n*B)
    expected <- matrix(x[donor] + h*sqrt(2/3)*qt(pnorm(noise), 6), n, B)
    end <- get(".Random.seed", .GlobalEnv)
    for (chunk in c(1L, 7L, 32L)) {
      tasks <- .np_udens_q6_plan(start, n, B, chunk)
      actual <- do.call(cbind, lapply(tasks, .np_udens_q6_samples, x=x, h=h))
      expect_identical(actual, expected)
      expect_identical(attr(tasks, "rng_final_state"), end)
      expect_identical(get(".Random.seed", .GlobalEnv), end)
    }
    reference <- .np_udens_q6_reference(x, h, c(-1, 0, 1), "uniform")
    full <- .np_udens_q6_maximum(expected, reference)$maximum
    parts <- lapply(tasks, .np_udens_q6_chunk, reference=reference)
    expect_equal(do.call(rbind, parts)[, 1L], full, tolerance=1e-12)
    tasks[[1L]]$noise.end[2L] <- tasks[[1L]]$noise.end[2L]+1L
    expect_error(.np_udens_q6_samples(tasks[[1L]], x, h), "inconsistent streams")
    expect_identical(get(".Random.seed", .GlobalEnv), end)
    expect_error(.np_udens_q6_plan(start, n, 0L), "invalid Q6 replicate")
  })
})

test_that("Q6 admission is an explicit isolated new law", {
  a <- list(errors="bootstrap", bootstrap="smoothed-q6", center="bias-corrected",
    band="simultaneous", boot.control=np_boot_control(calibration="symmetric-sd"), B=199L)
  expect_true(.np_udens_q6_options(a, TRUE))
  expect_error(.np_udens_q6_options(a), "only.*density")
  for (key in c("errors", "center", "band", "boot.control")) {
    b <- a; b[[key]] <- NULL
    expect_error(.np_udens_q6_options(b, TRUE), "requires explicit")
  }
  for (band in c("pmzsd")) {
    b<-a; b$band<-band; expect_error(.np_udens_q6_options(b,TRUE), "requires explicit")
  }
  for (cal in c("default", "absolute", "smoothed-percentile")) {
    b<-a; b$boot.control<-np_boot_control(calibration=cal)
    expect_error(.np_udens_q6_options(b,TRUE), "requires explicit")
  }
  for (B in c(.5, 0, Inf, NA_real_, .Machine$integer.max+1)) {
    b<-a;b$B<-B;expect_error(.np_udens_q6_options(b,TRUE), "representable integer")
  }
  for (extra in list(list(bias.control=np_udens_bias_control()),list(gradients=TRUE),
      list(gaussian.draws=199L),list(boot.control=np_boot_control(blocklen=2L,calibration="symmetric-sd")),
      list(boot.control=np_boot_control(nonfixed="frozen",calibration="symmetric-sd")),
      list(boot.control=np_boot_control(wild="rademacher",calibration="symmetric-sd"))))
    expect_error(.np_udens_q6_options(modifyList(a,extra),TRUE), "does not accept")
  b <- structure(list(type="fixed",ckerbound="none",ckertype="uniform",
    ckerorder=2L,ndim=1L,ncon=1L),class="bandwidth")
  expect_null(.np_udens_q6_spec(b))
  for (change in list(list(type="generalized_nn"),list(type="adaptive_nn"),
      list(ckerbound="fixed"),list(ckerbound="range"),list(ckerorder=4L),
      list(ndim=0L,ncon=0L),list(ndim=2L,ncon=2L),list(ndim=3L,ncon=3L),
      list(ncon=0L),list(ckertype="gaussian"),list(ckertype="beta")))
    expect_error(.np_udens_q6_spec(modifyList(b,change)), "requires q=1")
  expect_error(.np_plot_normalize_public_dots(a), "only.*density")
  norm <- .np_plot_normalize_public_dots(a,allow_smoothed_q6=TRUE)
  expect_identical(norm$plot.errors.boot.method,"smoothed-q6")
  expect_identical(norm$plot.errors.boot.calibration,"symmetric-sd")
  expect_identical(.np_plot_normalize_public_dots(norm,allow_smoothed_q6=TRUE),norm)
  # Adding this enum cannot admit it to another common-option caller.
  common <- list(plot.behavior="data",plot.errors.method="bootstrap",plot.errors.boot.method="smoothed-q6",
    plot.errors.boot.blocklen=NULL,plot.errors.center="bias-corrected",plot.errors.type="simultaneous",
    plot.errors.alpha=.05,plot.errors.style="band",plot.errors.bar="|",xdat=data.frame(x=1:10),common.scale=TRUE,ylim=NULL)
  expect_error(do.call(.np_plot_normalize_common_options,common), "arg")
  expect_identical(do.call(.np_plot_normalize_common_options,c(common,list(allow_smoothed_q6=TRUE)))$plot.errors.boot.method,"smoothed-q6")
  expect_false(.np_udens_q6_active(NULL))
})


test_that("Q6 scalar normalization cannot bypass its explicit admission", {
  a <- list(errors="bootstrap", center="bias-corrected", band="simultaneous",
    boot.control=np_boot_control(calibration="symmetric-sd"), B=199L)
  for (law in list(factor("smoothed-q6"), c(named="smoothed-q6"))) {
    a$bootstrap <- law
    z <- .np_plot_normalize_public_dots(a, allow_smoothed_q6=TRUE)
    expect_identical(z$plot.errors.boot.method, "smoothed-q6")
    expect_error(.np_plot_normalize_public_dots(a), "only.*density")
    for (key in c("errors", "center", "band", "boot.control")) {
      b <- a; b[[key]] <- NULL
      expect_error(.np_plot_normalize_public_dots(b, allow_smoothed_q6=TRUE), "requires explicit")
    }
    b <- a; b$center <- "estimate"
    expect_error(.np_plot_normalize_public_dots(b, allow_smoothed_q6=TRUE), "requires explicit")
    b <- a; b$B <- 1.5
    expect_error(.np_plot_normalize_public_dots(b, allow_smoothed_q6=TRUE), "representable integer")
    b <- a; b$gaussian.draws <- 199L
    expect_error(.np_plot_normalize_public_dots(b, allow_smoothed_q6=TRUE), "does not accept")
  }
})
