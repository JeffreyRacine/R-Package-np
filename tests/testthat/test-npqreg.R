test_that("npqreg basic functionality works", {
  data("cps71")
  cps71_sub <- cps71[1:50, ]
  
  # Quantile regression needs a condbandwidth object
  bw <- npcdistbw(xdat=cps71_sub$age, ydat=cps71_sub$logwage, 
                  bws=c(0.5, 5), bandwidth.compute=FALSE)
  
  # Median regression
  model <- npqreg(bws=bw, tau=0.5)
  
  expect_s3_class(model, "qregression")
  expect_type(predict(model), "double")
  expect_equal(length(predict(model)), 50)
  expect_equal(model$tau, 0.5)
  expect_error(gradients(model), "fit the model with gradients=TRUE", fixed = TRUE)
  expect_error(gradients(model, errors = TRUE), "fit the model with gradients=TRUE", fixed = TRUE)
  
  expect_output(summary(model))
})

test_that("npqreg works with multiple taus", {
  data("cps71")
  cps71_sub <- cps71[1:30, ]
  bw <- npcdistbw(xdat=cps71_sub$age, ydat=cps71_sub$logwage, 
                  bws=c(0.5, 5), bandwidth.compute=FALSE)
  
  # npqreg only takes a single tau at a time according to some versions, 
  # but let's see if it works as extra arg or in bws.
  # Actually, npqreg usage says tau is an argument.
  
  model_q25 <- npqreg(bws=bw, tau=0.25)
  expect_equal(model_q25$tau, 0.25)
})

test_that("npqreg formula route strips fit controls before bandwidth selection", {
  set.seed(20260511)
  d <- data.frame(x = runif(36), y = rnorm(36))

  fit <- npqreg(y ~ x, data = d, tau = c(0.25, 0.5),
                gradients = TRUE, tol = 1e-4, small = 1e-5,
                itmax = 1000L, nmulti = 1L)

  expect_s3_class(fit, "qregression")
  expect_equal(fit$tau, c(0.25, 0.5))
  expect_true(fit$gradients)
  expect_equal(ncol(fitted(fit)), 2L)
})

test_that("npqreg formula route honors manual bws with gradients", {
  set.seed(20260513)
  d <- data.frame(x = runif(35), y = rnorm(35))

  fit <- npqreg(y ~ x, data = d, tau = c(0.25, 0.5),
                gradients = TRUE, bws = c(0.3, 0.3),
                bandwidth.compute = FALSE)

  expect_s3_class(fit, "qregression")
  expect_equal(dim(fitted(fit)), c(nrow(d), 2L))
  expect_equal(dim(gradients(fit)), c(nrow(d), 1L, 2L))
})

test_that("npqreg validates tau and continuous responses before fitting", {
  set.seed(20260513)
  d <- data.frame(x = runif(24), y = rnorm(24))
  bw <- npcdistbw(
    xdat = d["x"],
    ydat = d["y"],
    bws = c(0.4, 0.4),
    bandwidth.compute = FALSE
  )

  expect_error(
    npqreg(bws = bw, txdat = d["x"], tydat = d$y, tau = 0),
    "'tau' must contain numeric values in (0,1)",
    fixed = TRUE
  )
  expect_error(
    npqreg(bws = bw, txdat = d["x"], tydat = d$y, tau = c(0.25, 1)),
    "'tau' must contain numeric values in (0,1)",
    fixed = TRUE
  )
  expect_error(
    npqreg(bws = bw, txdat = d["x"], tydat = d$y, tau = NA_real_),
    "'tau' must contain numeric values in (0,1)",
    fixed = TRUE
  )
  expect_error(
    npqreg(bws = bw, txdat = d["x"], tydat = d$y, tau = "0.5"),
    "'tau' must contain numeric values in (0,1)",
    fixed = TRUE
  )

  d$y_factor <- factor(rep(0:1, length.out = nrow(d)))
  expect_error(
    npqreg(y_factor ~ x, data = d, bws = c(0.4, 0.4),
           bandwidth.compute = FALSE),
    "'tydat' is not continuous",
    fixed = TRUE
  )
})

test_that("npqreg vector-tau plots support legend controls", {
  set.seed(20260511)
  n <- 24L
  xdat <- data.frame(x = seq(0.05, 0.95, length.out = n))
  ydat <- sin(2 * pi * xdat$x) + rnorm(n, sd = 0.05)
  bw <- npcdistbw(
    xdat = xdat,
    ydat = ydat,
    bws = c(0.5, 0.5),
    bandwidth.compute = FALSE
  )
  fit <- npqreg(bws = bw, tau = c(0.25, 0.5, 0.75))

  old.dev <- grDevices::dev.cur()
  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit({
    grDevices::dev.off()
    if (old.dev > 1L)
      grDevices::dev.set(old.dev)
  }, add = TRUE)

  expect_silent(plot(fit, perspective = FALSE, errors = "none", neval = 6L,
                     legend = NULL))
  expect_silent(plot(fit, perspective = FALSE, errors = "none", neval = 6L,
                     legend = NA))
  expect_silent(plot(fit, perspective = FALSE, errors = "none", neval = 6L,
                     legend = "bottomright"))
  expect_silent(plot(fit, perspective = FALSE, errors = "none", neval = 6L,
                     legend = list(x = "bottomleft",
                                   legend = c("low", "mid", "high"),
                                   bty = "o")))
  expect_error(
    plot(fit, perspective = FALSE, errors = "none", neval = 6L,
         legend = 1),
    "legend must be TRUE/FALSE"
  )
})

test_that("npqreg plot preserves NOMAD-selected LP degree metadata", {
  set.seed(70511)
  n <- 24L
  d <- data.frame(
    x = runif(n),
    y = sin(runif(n)) + rnorm(n, sd = 0.05)
  )
  bw <- npcdistbw(
    y ~ x,
    data = d,
    nomad = TRUE,
    search.engine = "nomad",
    degree.min = 0L,
    degree.max = 1L,
    degree.start = 1L,
    degree.verify = FALSE,
    nmulti = 1L,
    nomad.nmulti = 0L
  )
  expect_identical(bw$regtype, "lp")
  expect_identical(bw$regtype.engine, "lp")
  expect_false(is.null(bw$degree))
  expect_false(is.null(bw$degree.engine))

  fit <- npqreg(bws = bw, tau = 0.5)
  out <- plot(
    fit,
    output = "data",
    perspective = FALSE,
    tau = c(0.25, 0.5),
    errors = "none",
    neval = 5L
  )

  expect_true(all(vapply(out, inherits, logical(1), "qregression")))
  expect_true(all(vapply(out, function(xi) identical(xi$bws$regtype.engine, "lp"), logical(1))))
  expect_true(all(vapply(out, function(xi) identical(xi$bws$degree.engine, bw$degree.engine), logical(1))))
})

test_that("npqreg gradients populate qregression objects", {
  set.seed(70507)
  n <- 150L
  xdat <- data.frame(
    x1 = runif(n, -1, 1),
    x2 = runif(n, -1, 1)
  )
  ydat <- xdat$x1 - 0.5 * xdat$x2 + rnorm(n, sd = 0.15)
  exdat <- data.frame(
    x1 = seq(-0.40, 0.40, length.out = 5L),
    x2 = seq(0.35, -0.35, length.out = 5L)
  )
  bw <- npcdistbw(
    xdat = xdat,
    ydat = ydat,
    bws = c(0.80, 0.80, 0.80),
    bandwidth.compute = FALSE
  )

  fit <- npqreg(bws = bw, txdat = xdat, tydat = ydat, exdat = exdat,
                tau = 0.45, gradients = TRUE, tol = 1e-6, small = 1e-7)
  ref <- getFromNamespace(".np_plot_quantile_eval", "np")(
    bws = bw,
    txdat = xdat,
    tydat = ydat,
    exdat = exdat,
    tau = 0.45,
    gradients = TRUE,
    tol = 1e-6,
    small = 1e-7
  )

  expect_s3_class(fit, "qregression")
  expect_true(fit$gradients)
  expect_equal(dim(fit$quantgrad), c(nrow(exdat), ncol(xdat)))
  expect_equal(dim(fit$quantgerr), c(nrow(exdat), ncol(xdat)))
  expect_true(all(is.finite(fit$quantgrad)))
  expect_true(all(is.finite(fit$quantgerr)))
  expect_true(all(is.finite(se(fit))))
  expect_true(all(se(fit) > 0))
  expect_equal(fit$quantgrad, ref$quantgrad, tolerance = 0)
  expect_equal(gradients(fit), fit$quantgrad, tolerance = 0)
  expect_equal(gradients(fit, errors = TRUE), fit$quantgerr, tolerance = 0)

  qdelta <- getFromNamespace(".npqreg_quantile_delta_from_conditional", "np")(
    bws = bw,
    xdat = xdat,
    ydat = ydat,
    exdat = exdat,
    quantile = fitted(fit),
    gradients = TRUE
  )
  expect_equal(se(fit), qdelta$quanterr, tolerance = 0)
  expect_equal(fit$quantgrad, qdelta$quantgrad, tolerance = 0)
  expect_equal(fit$quantgerr, qdelta$quantgerr, tolerance = 0)

  eps <- 1e-2
  for (j in seq_along(exdat)) {
    xp <- xm <- exdat
    xp[[j]] <- xp[[j]] + eps
    xm[[j]] <- xm[[j]] - eps
    qp <- fitted(npqreg(bws = bw, txdat = xdat, tydat = ydat, exdat = xp,
                        tau = 0.45, gradients = FALSE, tol = 1e-6, small = 1e-7))
    qm <- fitted(npqreg(bws = bw, txdat = xdat, tydat = ydat, exdat = xm,
                        tau = 0.45, gradients = FALSE, tol = 1e-6, small = 1e-7))
    fd <- (qp - qm) / (2 * eps)
    expect_lt(max(abs(fit$quantgrad[, j] - fd), na.rm = TRUE), 1e-3)
  }
})

test_that("npqreg gradient plots support asymptotic errors", {
  set.seed(70512)
  n <- 48L
  xdat <- data.frame(x = sort(runif(n)))
  ydat <- sin(2 * pi * xdat$x) + rnorm(n, sd = 0.15)
  bw <- npcdistbw(
    xdat = xdat,
    ydat = ydat,
    bws = c(0.45, 0.45),
    bandwidth.compute = FALSE,
    cxkerbound = "range",
    cykerbound = "range"
  )
  fit <- npqreg(bws = bw, tau = c(0.25, 0.5, 0.75), gradients = TRUE)

  out <- plot(
    fit,
    gradients = TRUE,
    output = "data",
    perspective = FALSE,
    errors = "asymptotic",
    neval = 8L
  )

  expect_true(all(vapply(out, inherits, logical(1), "qregression")))
  expect_true(all(is.finite(gradients(fit))))
  expect_true(all(is.finite(gradients(fit, errors = TRUE))))
})
