test_that("partial-linear covariance work is controlled at its shared owner", {
  ns <- asNamespace(getNamespaceName(environment(npplreg)))
  solve.linear <- get(".np_plreg_linear_solve", ns)
  X <- cbind(x = c(-2, -1, 1, 2, -2, -1, 1, 2),
             w = c(-1, 1, -1, 1, 1, -1, 1, -1))
  y <- as.vector(X %*% c(2, -.3)) + c(.1, -.2, .2, -.1, -.1, .2, -.2, .1)
  on <- solve.linear(y, X, y, rep(0, length(y)), 1L, se = TRUE)
  off <- solve.linear(y, X, y, rep(0, length(y)), 1L, se = FALSE)
  expect_identical(off$coef, on$coef)
  expect_identical(off$train.fit, on$train.fit)
  expect_identical(off$qr, on$qr)
  expect_null(off$vcov)
  expect_null(off$se)

  # Instrument a private closure, not the package's production namespace.
  guarded <- solve.linear
  environment(guarded) <- new.env(parent = environment(solve.linear))
  environment(guarded)$chol2inv <- function(...)
    stop("unexpected covariance inverse", call. = FALSE)
  expect_identical(guarded(y, X, y, rep(0, length(y)), 1L, se = FALSE), off)
  expect_error(guarded(y, X, y, rep(0, length(y)), 1L, se = TRUE),
               "unexpected covariance inverse", fixed = TRUE)
  expect_error(guarded(y, cbind(X[, 1L], X[, 1L]), y, rep(0, length(y)),
                       1L, se = FALSE), "rank deficient", fixed = TRUE)
})

test_that("partial-linear uncertainty extraction is explicit and historical-state aware", {
  historical <- structure(list(xcoef = c(x = 2), xcoeferr = c(x = NA_real_),
                               xcoefvcov = matrix(0, 1L, 1L)),
                          class = "plregression")
  expect_identical(coef(historical, se = TRUE), historical$xcoeferr)
  expect_identical(vcov(historical), historical$xcoefvcov)
  omitted <- historical
  omitted$se <- FALSE
  expect_error(coef(omitted, se = TRUE),
               "npplreg(bws = omitted$bws, se = TRUE)", fixed = TRUE)
  expect_error(vcov(omitted), "without repeating bandwidth search", fixed = TRUE)
  omitted$bw <- structure(list(), class = "plbandwidth")
  expect_error(vcov(omitted), "npplreg(bws = omitted$bw, se = TRUE)", fixed = TRUE)
  expect_identical(coef(omitted), historical$xcoef)

  for (method in c("formula", "default", "plbandwidth"))
    expect_identical(formals(getS3method("npplreg", method))$se, FALSE)
  for (method in c("formula", "default", "plbandwidth")) {
    dispatch <- getS3method("npplreg", method)
    before.dots <- match("...", names(formals(dispatch))) - 1L
    supplied <- c(list(as.name("npplreg")), rep(list(NULL), before.dots),
                  list(quote(legacy.positional.argument)))
    matched <- match.call(dispatch, as.call(supplied), expand.dots = FALSE)
    expect_identical(matched[["..."]][[1L]], quote(legacy.positional.argument))
    expect_null(matched[["se"]])
  }
  expect_error(npplreg(se = NA), "'se' must", fixed = TRUE)
  expect_error(npplreg(s.e = stop("should not be evaluated")),
               "did you mean 'se'", fixed = TRUE)
})

test_that("partial-linear public fitting keeps points and optionally stores covariance", {
  if (exists("spawn_mpi_slaves", mode = "function")) {
    if (!spawn_mpi_slaves()) skip("Could not initialize MPI context")
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }
  dat <- expand.grid(z = letters[1:2], x = c(-2, -1, 1, 2), e = c(-1, 1))
  dat$z <- factor(dat$z)
  dat$y <- 2 * dat$x + as.numeric(dat$z) + .1 * dat$e
  bw <- npplregbw(y ~ x | z, data = dat, bws = matrix(0, 2L, 1L),
                  bandwidth.compute = FALSE)
  plain <- npplreg(bws = bw)
  full <- npplreg(bws = bw, se = TRUE)
  expect_identical(plain$se, FALSE)
  expect_identical(full$se, TRUE)
  expect_null(plain$xcoeferr)
  expect_null(plain$xcoefvcov)
  expect_identical(coef(plain), coef(full))
  expect_identical(fitted(plain), fitted(full))
  expect_identical(coef(full, se = TRUE), sqrt(diag(vcov(full))))
  expect_error(vcov(plain), "npplreg(bws = plain$bws, se = TRUE)", fixed = TRUE)
})
