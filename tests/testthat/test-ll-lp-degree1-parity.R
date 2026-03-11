test_that("npreg cv objective and bandwidths match for ll and lp(degree=1)", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260227)
  x <- runif(70)
  y <- sin(2 * pi * x) + rnorm(70, sd = 0.05)
  tx <- data.frame(x = x)

  for (m in c("cv.ls", "cv.aic")) {
    set.seed(90210)
    bw.ll <- npregbw(
      xdat = tx,
      ydat = y,
      regtype = "ll",
      bwmethod = m,
      nmulti = 1L
    )
    set.seed(90210)
    bw.lp <- npregbw(
      xdat = tx,
      ydat = y,
      regtype = "lp",
      basis = "glp",
      degree = 1L,
      bwmethod = m,
      nmulti = 1L
    )

    expect_equal(as.numeric(bw.ll$fval), as.numeric(bw.lp$fval), tolerance = 1e-10)
    expect_equal(as.numeric(bw.ll$bw), as.numeric(bw.lp$bw), tolerance = 1e-9)
  }
})

test_that("npreg and npreghat match for ll and lp(degree=1) in 1D", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260227)
  n <- 140
  x <- sort(runif(n))
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.02)
  tx <- data.frame(x = x)
  ex <- data.frame(x = seq(min(x), max(x), length.out = 50))

  bw.ll <- npregbw(
    xdat = tx,
    ydat = y,
    bws = 0.2,
    regtype = "ll",
    bandwidth.compute = FALSE
  )
  bw.lp <- npregbw(
    xdat = tx,
    ydat = y,
    bws = 0.2,
    regtype = "lp",
    basis = "glp",
    degree = 1L,
    bandwidth.compute = FALSE
  )

  fit.ll <- npreg(
    txdat = tx,
    tydat = y,
    exdat = ex,
    bws = bw.ll,
    gradients = TRUE,
    gradient.order = 1L,
    warn.glp.gradient = FALSE
  )
  fit.lp <- npreg(
    txdat = tx,
    tydat = y,
    exdat = ex,
    bws = bw.lp,
    gradients = TRUE,
    gradient.order = 1L,
    warn.glp.gradient = FALSE
  )

  expect_equal(as.numeric(fit.ll$mean), as.numeric(fit.lp$mean), tolerance = 1e-10)
  expect_equal(as.numeric(fit.ll$grad), as.numeric(fit.lp$grad), tolerance = 1e-10)

  H0.ll <- npreghat(bws = bw.ll, txdat = tx, exdat = ex)
  H0.lp <- npreghat(bws = bw.lp, txdat = tx, exdat = ex)
  H1.ll <- npreghat(bws = bw.ll, txdat = tx, exdat = ex, s = 1L)
  H1.lp <- npreghat(bws = bw.lp, txdat = tx, exdat = ex, s = 1L)

  expect_equal(as.numeric(H0.ll %*% y), as.numeric(H0.lp %*% y), tolerance = 1e-10)
  expect_equal(as.numeric(H1.ll %*% y), as.numeric(H1.lp %*% y), tolerance = 1e-10)
})

test_that("npscoef cv and estimation match for ll and lp(degree=1) in 1D", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260227)
  n <- 80
  x <- runif(n)
  z <- runif(n)
  y <- (0.5 + sin(2 * pi * x)) * z + rnorm(n, sd = 0.04)

  set.seed(7701)
  bw.ll.cv <- npscoefbw(
    xdat = x,
    zdat = z,
    ydat = y,
    regtype = "ll",
    nmulti = 1L,
    cv.iterate = FALSE
  )
  set.seed(7701)
  bw.lp.cv <- npscoefbw(
    xdat = x,
    zdat = z,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = 1L,
    nmulti = 1L,
    cv.iterate = FALSE
  )

  expect_equal(as.numeric(bw.ll.cv$fval), as.numeric(bw.lp.cv$fval), tolerance = 1e-10)
  expect_equal(as.numeric(bw.ll.cv$bw), as.numeric(bw.lp.cv$bw), tolerance = 1e-9)

  bw.ll <- npscoefbw(
    xdat = x,
    zdat = z,
    ydat = y,
    regtype = "ll",
    bws = 0.2,
    bandwidth.compute = FALSE
  )
  bw.lp <- npscoefbw(
    xdat = x,
    zdat = z,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = 1L,
    bws = 0.2,
    bandwidth.compute = FALSE
  )

  fit.ll <- npscoef(
    bws = bw.ll,
    txdat = data.frame(x = x),
    tzdat = data.frame(z = z),
    tydat = y,
    errors = FALSE
  )
  fit.lp <- npscoef(
    bws = bw.lp,
    txdat = data.frame(x = x),
    tzdat = data.frame(z = z),
    tydat = y,
    errors = FALSE
  )

  expect_equal(as.numeric(fit.ll$mean), as.numeric(fit.lp$mean), tolerance = 1e-10)
})

test_that("npindex cv, fit, and npindexhat(s=1) match for ll and lp(degree=1) in 1D", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260227)
  n <- 100
  x <- runif(n)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.04)
  tx <- data.frame(x = x)

  set.seed(8801)
  bw.ll.cv <- suppressWarnings(npindexbw(
    xdat = tx,
    ydat = y,
    regtype = "ll",
    method = "ichimura",
    nmulti = 1L
  ))
  set.seed(8801)
  bw.lp.cv <- suppressWarnings(npindexbw(
    xdat = tx,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = 1L,
    method = "ichimura",
    nmulti = 1L
  ))

  expect_equal(as.numeric(bw.ll.cv$fval), as.numeric(bw.lp.cv$fval), tolerance = 1e-10)
  expect_equal(as.numeric(bw.ll.cv$bw), as.numeric(bw.lp.cv$bw), tolerance = 1e-9)

  bw.ll <- suppressWarnings(npindexbw(
    xdat = tx,
    ydat = y,
    bws = c(1, 0.25),
    regtype = "ll",
    method = "ichimura",
    bandwidth.compute = FALSE
  ))
  bw.lp <- suppressWarnings(npindexbw(
    xdat = tx,
    ydat = y,
    bws = c(1, 0.25),
    regtype = "lp",
    basis = "glp",
    degree = 1L,
    method = "ichimura",
    bandwidth.compute = FALSE
  ))

  fit.ll <- npindex(bws = bw.ll, txdat = tx, tydat = y, exdat = data.frame(x = seq(0, 1, length.out = 40)))
  fit.lp <- npindex(bws = bw.lp, txdat = tx, tydat = y, exdat = data.frame(x = seq(0, 1, length.out = 40)))
  expect_equal(as.numeric(fit.ll$mean), as.numeric(fit.lp$mean), tolerance = 1e-10)

  H0.ll <- npindexhat(bws = bw.ll, txdat = tx)
  H0.lp <- npindexhat(bws = bw.lp, txdat = tx)
  H1.ll <- npindexhat(bws = bw.ll, txdat = tx, s = 1L)
  H1.lp <- npindexhat(bws = bw.lp, txdat = tx, s = 1L)
  expect_equal(as.numeric(H0.ll %*% y), as.numeric(H0.lp %*% y), tolerance = 1e-10)
  expect_equal(as.numeric(H1.ll %*% y), as.numeric(H1.lp %*% y), tolerance = 1e-10)
})

test_that("npplreg cv and estimation match for ll and lp(degree=1) in 1D", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260227)
  n <- 120
  z <- runif(n)
  x <- runif(n)
  y <- sin(2 * pi * z) + 1.2 * x + rnorm(n, sd = 0.04)

  set.seed(9901)
  bw.ll.cv <- npplregbw(
    xdat = data.frame(x = x),
    zdat = data.frame(z = z),
    ydat = y,
    regtype = "ll",
    nmulti = 1L
  )
  set.seed(9901)
  bw.lp.cv <- npplregbw(
    xdat = data.frame(x = x),
    zdat = data.frame(z = z),
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = 1L,
    nmulti = 1L
  )

  expect_equal(as.numeric(bw.ll.cv$fval), as.numeric(bw.lp.cv$fval), tolerance = 1e-10)
  expect_equal(as.numeric(bw.ll.cv$bw$yzbw$bw), as.numeric(bw.lp.cv$bw$yzbw$bw), tolerance = 1e-9)
  expect_equal(as.numeric(bw.ll.cv$bw[[2]]$bw), as.numeric(bw.lp.cv$bw[[2]]$bw), tolerance = 1e-9)

  bw.ll <- npplregbw(
    xdat = data.frame(x = x),
    zdat = data.frame(z = z),
    ydat = y,
    regtype = "ll",
    bws = matrix(c(0.2, 0.2), nrow = 2),
    bandwidth.compute = FALSE
  )
  bw.lp <- npplregbw(
    xdat = data.frame(x = x),
    zdat = data.frame(z = z),
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = 1L,
    bws = matrix(c(0.2, 0.2), nrow = 2),
    bandwidth.compute = FALSE
  )

  ex <- data.frame(x = runif(30))
  ez <- data.frame(z = runif(30))
  fit.ll <- npplreg(
    bws = bw.ll,
    txdat = data.frame(x = x),
    tzdat = data.frame(z = z),
    tydat = y,
    exdat = ex,
    ezdat = ez
  )
  fit.lp <- npplreg(
    bws = bw.lp,
    txdat = data.frame(x = x),
    tzdat = data.frame(z = z),
    tydat = y,
    exdat = ex,
    ezdat = ez
  )

  expect_equal(as.numeric(fit.ll$mean), as.numeric(fit.lp$mean), tolerance = 1e-10)
  expect_equal(as.numeric(fit.ll$xcoef), as.numeric(fit.lp$xcoef), tolerance = 1e-10)
})

test_that("npscoef multivariate cv and estimation match for ll and lp(degree=1)", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260306)
  n <- 65
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  z1 <- runif(n, -1, 1)
  z2 <- runif(n, -1, 1)
  b1 <- 1 + 0.6 * z1 - 0.2 * z2
  b2 <- -0.5 + 0.4 * z1 + 0.3 * z2
  y <- b1 * x1 + b2 * x2 + rnorm(n, sd = 0.08)

  tx <- data.frame(x1 = x1, x2 = x2)
  tz <- data.frame(z1 = z1, z2 = z2)

  set.seed(7702)
  bw.ll.cv <- npscoefbw(
    xdat = tx,
    zdat = tz,
    ydat = y,
    regtype = "ll",
    nmulti = 1L,
    cv.iterate = FALSE,
    backfit.iterate = FALSE
  )
  set.seed(7702)
  bw.lp.cv <- npscoefbw(
    xdat = tx,
    zdat = tz,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = c(1L, 1L),
    nmulti = 1L,
    cv.iterate = FALSE,
    backfit.iterate = FALSE
  )

  expect_equal(as.numeric(bw.ll.cv$fval), as.numeric(bw.lp.cv$fval), tolerance = 1e-10)
  expect_equal(as.numeric(bw.ll.cv$bw), as.numeric(bw.lp.cv$bw), tolerance = 1e-9)
  expect_equal(as.numeric(bw.ll.cv$num.feval), as.numeric(bw.lp.cv$num.feval), tolerance = 0)

  bw.ll <- npscoefbw(
    xdat = tx,
    zdat = tz,
    ydat = y,
    regtype = "ll",
    bws = c(0.7, 0.7),
    bandwidth.compute = FALSE
  )
  bw.lp <- npscoefbw(
    xdat = tx,
    zdat = tz,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = c(1L, 1L),
    bws = c(0.7, 0.7),
    bandwidth.compute = FALSE
  )

  fit.ll <- npscoef(bws = bw.ll, txdat = tx, tzdat = tz, tydat = y, errors = FALSE)
  fit.lp <- npscoef(bws = bw.lp, txdat = tx, tzdat = tz, tydat = y, errors = FALSE)

  expect_equal(as.numeric(fit.ll$mean), as.numeric(fit.lp$mean), tolerance = 1e-10)
  expect_equal(as.numeric(fit.ll$beta), as.numeric(fit.lp$beta), tolerance = 1e-10)
})

test_that("npindex multivariate cv and estimation match for ll and lp(degree=1)", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  n <- 75

  run_case <- function(method) {
    set.seed(if (identical(method, "ichimura")) 20260307 else 20260308)
    x1 <- runif(n, -1, 1)
    x2 <- runif(n, -1, 1)
    x3 <- runif(n, -1, 1)
    eta <- x1 + 0.5 * x2 - 0.4 * x3
    y <- if (identical(method, "ichimura")) {
      eta + 0.3 * eta^2 + rnorm(n, sd = 0.1)
    } else {
      rbinom(n, size = 1L, prob = plogis(eta))
    }
    tx <- data.frame(x1 = x1, x2 = x2, x3 = x3)

    set.seed(8802)
    bw.ll.cv <- suppressWarnings(npindexbw(
      xdat = tx,
      ydat = y,
      regtype = "ll",
      method = method,
      nmulti = 1L
    ))
    set.seed(8802)
    bw.lp.cv <- suppressWarnings(npindexbw(
      xdat = tx,
      ydat = y,
      regtype = "lp",
      basis = "glp",
      degree = 1L,
      method = method,
      nmulti = 1L
    ))

    expect_equal(as.numeric(bw.ll.cv$fval), as.numeric(bw.lp.cv$fval), tolerance = 1e-10)
    expect_equal(
      c(as.numeric(bw.ll.cv$beta), as.numeric(bw.ll.cv$bw)),
      c(as.numeric(bw.lp.cv$beta), as.numeric(bw.lp.cv$bw)),
      tolerance = 1e-9
    )
    expect_equal(as.numeric(bw.ll.cv$num.feval), as.numeric(bw.lp.cv$num.feval), tolerance = 0)

    bws.fixed <- c(1, 0.3, -0.2, 0.6)
    bw.ll <- suppressWarnings(npindexbw(
      xdat = tx,
      ydat = y,
      bws = bws.fixed,
      regtype = "ll",
      method = method,
      bandwidth.compute = FALSE
    ))
    bw.lp <- suppressWarnings(npindexbw(
      xdat = tx,
      ydat = y,
      bws = bws.fixed,
      regtype = "lp",
      basis = "glp",
      degree = 1L,
      method = method,
      bandwidth.compute = FALSE
    ))

    fit.ll <- npindex(bws = bw.ll, txdat = tx, tydat = y)
    fit.lp <- npindex(bws = bw.lp, txdat = tx, tydat = y)
    expect_equal(as.numeric(fit.ll$mean), as.numeric(fit.lp$mean), tolerance = 1e-10)
    expect_equal(as.numeric(fit.ll$beta), as.numeric(fit.lp$beta), tolerance = 1e-10)
  }

  run_case("ichimura")
  run_case("kleinspady")
})

test_that("npplreg multivariate cv and estimation match for ll and lp(degree=1)", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260309)
  n <- 70
  z1 <- runif(n, -1, 1)
  z2 <- runif(n, -1, 1)
  x1 <- 0.8 * z1 + rnorm(n, sd = 0.2)
  x2 <- -0.6 * z2 + rnorm(n, sd = 0.2)
  y <- 1 + 1.2 * x1 - 0.7 * x2 + 0.4 * z1^2 - 0.3 * z1 * z2 + rnorm(n, sd = 0.1)

  tx <- data.frame(x1 = x1, x2 = x2)
  tz <- data.frame(z1 = z1, z2 = z2)

  set.seed(9902)
  bw.ll.cv <- npplregbw(
    xdat = tx,
    zdat = tz,
    ydat = y,
    regtype = "ll",
    nmulti = 1L
  )
  set.seed(9902)
  bw.lp.cv <- npplregbw(
    xdat = tx,
    zdat = tz,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = c(1L, 1L),
    nmulti = 1L
  )

  expect_equal(as.numeric(bw.ll.cv$fval), as.numeric(bw.lp.cv$fval), tolerance = 1e-10)
  expect_equal(as.numeric(bw.ll.cv$bw$yzbw$bw), as.numeric(bw.lp.cv$bw$yzbw$bw), tolerance = 1e-9)
  expect_equal(as.numeric(bw.ll.cv$bw[[2L]]$bw), as.numeric(bw.lp.cv$bw[[2L]]$bw), tolerance = 1e-9)
  expect_equal(as.numeric(bw.ll.cv$bw[[3L]]$bw), as.numeric(bw.lp.cv$bw[[3L]]$bw), tolerance = 1e-9)
  expect_equal(as.numeric(bw.ll.cv$num.feval), as.numeric(bw.lp.cv$num.feval), tolerance = 0)

  bws.fixed <- matrix(c(0.3, 0.3, 0.3, 0.3, 0.3, 0.3), nrow = 3L)
  bw.ll <- npplregbw(
    xdat = tx,
    zdat = tz,
    ydat = y,
    regtype = "ll",
    bws = bws.fixed,
    bandwidth.compute = FALSE
  )
  bw.lp <- npplregbw(
    xdat = tx,
    zdat = tz,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = c(1L, 1L),
    bws = bws.fixed,
    bandwidth.compute = FALSE
  )

  fit.ll <- npplreg(bws = bw.ll, txdat = tx, tzdat = tz, tydat = y)
  fit.lp <- npplreg(bws = bw.lp, txdat = tx, tzdat = tz, tydat = y)

  expect_equal(as.numeric(fit.ll$mean), as.numeric(fit.lp$mean), tolerance = 1e-10)
  expect_equal(as.numeric(fit.ll$xcoef), as.numeric(fit.lp$xcoef), tolerance = 1e-10)
})
