test_that("generic adaptive LP CV fallback retains delete-one geometry", {
  skip_on_cran()
  if (!spawn_mpi_slaves(1L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE, np.tree = FALSE, np.extendednn = FALSE)
  on.exit(options(old), add = TRUE)
  injection <- "NP_RMPI_INJECT_REG_ADAPTIVE_BLAS_FAIL_RANK"
  old.injection <- Sys.getenv(injection, unset = NA_character_)
  on.exit(if (is.na(old.injection)) Sys.unsetenv(injection) else
    do.call(Sys.setenv, setNames(list(old.injection), injection)), add = TRUE)
  x <- as.data.frame(lapply(1:2, function(j) c(rep(0, 32L),
    sin((1:16)*sqrt(j+1)) + (1:16)/(16*(j+1)))))
  y <- sin((1:48)/3) + (1:48)/90
  raw <- function(k, degree) {
    b <- npregbw(xdat = x, ydat = y, bws = k, bwtype = "adaptive_nn",
      regtype = "lp", degree = degree, nomad = FALSE,
      bandwidth.compute = FALSE, ckertype = "epanechnikov")
    .npregbw_eval_only(x, y, b, invalid.penalty = "dbmax")$objective
  }
  for (inject in c("", "all")) {
    if (nzchar(inject)) do.call(Sys.setenv, setNames(list(inject), injection)) else
      Sys.unsetenv(injection)
    # The full-sample radius cannot replace a rejected delete-one radius.
    expect_identical(raw(c(6, 10), c(1L, 0L)), .Machine$double.xmax)
    expect_identical(raw(c(7, 7), c(1L, 1L)), .Machine$double.xmax)
  }
  # Reuse after rejection: a valid fallback must remain available, not be
  # turned into an unconditional failure or an invalid finite penalty.
  generic <- raw(c(40, 42), c(1L, 0L))
  Sys.unsetenv(injection)
  specialized <- raw(c(40, 42), c(1L, 0L))
  expect_true(is.finite(generic) && abs(generic) < .Machine$double.xmax)
  expect_equal(generic, specialized, tolerance = 2e-12)
  expect_true(is.finite(raw(c(40, 42), c(0L, 0L))))
})
