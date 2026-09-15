test_that("copula coordinates follow evaluation state, not published column names", {
  if (exists("spawn_mpi_slaves", mode = "function")) {
    if (!spawn_mpi_slaves()) skip("Could not initialize MPI context")
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }
  ns <- asNamespace(getNamespaceName(environment(npcopula)))
  coordinates <- get(".npcopula_eval_xgrid", ns)
  set.seed(71)
  plain <- data.frame(x = rnorm(25), y = rnorm(25))
  renamed <- list(c("x", "y"), c("x y", "x-y"), c("u1", "copula"))
  for (target in c("distribution", "density")) {
    make.bw <- if (target == "density") npudensbw else npudistbw
    for (grid in c(FALSE, TRUE)) {
      reference <- NULL
      for (nms in renamed) {
        dat <- plain
        names(dat) <- nms
        bw <- make.bw(dat = dat, bws = c(.8, .9), bandwidth.compute = FALSE)
        fit <- npcopula(bws = bw, data = dat, se = TRUE,
          u = if (grid) matrix(rep(seq(.1, .9, length.out = 5), 2), 5) else NULL,
          n.quasi.inv = 30L)
        retained <- fit
        xy <- coordinates(fit)
        expect_identical(names(xy), nms)
        # The independent oracle is the producer's physical block, not training
        # data merely because this grid also has 25 rows.
        oracle <- if (grid) as.data.frame(fit)[, 4:5, drop = FALSE] else dat
        expect_identical(unname(as.matrix(xy)), unname(as.matrix(oracle)))
        result <- list(point = fitted(fit), uncertainty = se(fit),
          geometry = unname(as.matrix(xy)))
        if (is.null(reference)) reference <- result else expect_identical(result, reference)
        expect_identical(fit, retained)
        # Serialize only the coordinate contract, not the test's calling
        # environment retained by bandwidth metadata.
        portable <- fit[c("evaluation", "xnames", "copula", "grid.dim", "data", "eval")]
        class(portable) <- "npcopula"
        expect_identical(coordinates(unserialize(serialize(portable, NULL))), xy)
        if (grid) {
          pd <- plot(fit, output = "data", errors = "asymptotic")
          expect_equal(nrow(pd), 25L)
          expect_identical(get(".npcopula_grid_eval", ns)(fit)$xgrid, xy)
        }
      }
    }
  }
})

test_that("transformed copula formula keeps grid SE and prediction geometry", {
  if (exists("spawn_mpi_slaves", mode = "function")) {
    if (!spawn_mpi_slaves()) skip("Could not initialize MPI context")
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }
  set.seed(42)
  for (n in c(25L, 28L)) {
    dat <- data.frame(x = runif(n), z = exp(rnorm(n)))
    plain <- data.frame(x = dat$x, logz = log(dat$z))
    named <- npcopula(~ x + log(z), data = dat, bws = c(.8, .8),
      bandwidth.compute = FALSE, neval = 5L, n.quasi.inv = 30L, se = TRUE)
    control <- npcopula(~ x + logz, data = plain, bws = c(.8, .8),
      bandwidth.compute = FALSE, neval = 5L, n.quasi.inv = 30L, se = TRUE)
    expect_identical(fitted(named), fitted(control))
    expect_identical(se(named), se(control))
    u <- data.frame(u1 = c(.2, .7), u2 = c(.3, .8))
    expect_identical(predict(named, newdata = u, se.fit = TRUE, n.quasi.inv = 30L),
      predict(control, newdata = u, se.fit = TRUE, n.quasi.inv = 30L))
    set.seed(100)
    a <- plot(named, output = "data", errors = "bootstrap", band = "pmzsd", B = 5L)
    rng <- .Random.seed
    set.seed(100)
    b <- plot(control, output = "data", errors = "bootstrap", band = "pmzsd", B = 5L)
    expect_identical(unname(as.matrix(a)), unname(as.matrix(b)))
    expect_identical(.Random.seed, rng)
  }
})

test_that("malformed copula grid metadata cannot select training coordinates", {
  ns <- asNamespace(getNamespaceName(environment(npcopula)))
  coordinates <- get(".npcopula_eval_xgrid", ns)
  # Deliberately equal row counts: no name- or row-count-based rescue.
  x <- structure(list(evaluation = "grid", xnames = c("x", "log(y)"),
    copula = c(.2, .4, .6, .8), grid.dim = c(2L, 2L),
    data = data.frame(x = 1:4, check.names = FALSE, "log(y)" = 5:8),
    eval = data.frame(copula = c(.2, .4, .6, .8), u1 = c(.2, .8, .2, .8),
      u2 = c(.2, .2, .8, .8), x = 9:12, log.y. = 13:16)),
    class = "npcopula")
  expect_identical(coordinates(x)[[1L]], 9:12)
  for (kind in c("missing-state", "missing-column", "bad-prefix", "bad-dim")) {
    bad <- x
    if (kind == "missing-state") bad$evaluation <- NULL
    if (kind == "missing-column") bad$eval <- bad$eval[, -5L]
    if (kind == "bad-prefix") names(bad$eval)[2L] <- "wrong"
    if (kind == "bad-dim") bad$grid.dim <- c(2L, 3L)
    expect_error(coordinates(bad), "evaluation coordinates.*refit", ignore.case = TRUE)
  }
})
