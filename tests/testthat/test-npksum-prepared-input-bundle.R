test_that("npksum bundles publications without changing values or cleanup", {
  skip_if_not(spawn_mpi_slaves())
  on.exit(close_mpi_slaves(), add = TRUE)
  old <- options(np.messages = FALSE, npRmpi.profile.level = "detailed",
                 npRmpi.autodispatch.arg.broadcast.threshold = 1L)
  on.exit(options(old), add = TRUE)
  set.seed(834L)
  dat <- data.frame(x = runif(60), f = factor(rep(letters[1:3], 20)), y = rnorm(60))
  fit <- npksum(y ~ x + f, data = dat, bws = c(.2, .3))
  direct <- npksum(txdat = dat[c("x", "f")], tydat = dat$y, bws = c(.2, .3))
  expect_identical(fit$ksum, direct$ksum)
  expect_false(any(grepl(".__npRmpi_autod_", deparse(fit$call), fixed = TRUE)))
  expect_identical(sum(fit$timing.profile$comm_notes == "mpi.bcast.Robj2slave"), 1L)
  expect_identical(sum(direct$timing.profile$comm_notes == "mpi.bcast.Robj2slave"), 1L)
  globals <- function() grep("^\\.__npRmpi_autod_", ls(.GlobalEnv, all.names = TRUE), value = TRUE)
  expect_length(globals(), 0L)
  remaining <- npRmpi:::mpi.remote.exec(
    grep("^\\.__npRmpi_autod_", ls(.GlobalEnv, all.names = TRUE), value = TRUE),
    simplify = FALSE)
  expect_true(all(lengths(remaining) == 0L))
  expect_error(npksum(txdat = dat[c("x", "f")], tydat = dat$y[-1], bws = c(.2, .3)),
               "number of explanatory data 'txdat' and dependent data 'tydat' do not match",
               fixed = TRUE)
  expect_length(globals(), 0L)
  remaining <- npRmpi:::mpi.remote.exec(
    grep("^\\.__npRmpi_autod_", ls(.GlobalEnv, all.names = TRUE), value = TRUE),
    simplify = FALSE)
  expect_true(all(lengths(remaining) == 0L))
  expect_identical(npksum(txdat = dat[c("x", "f")], tydat = dat$y,
                         bws = c(.2, .3))$ksum, direct$ksum)

  options(npRmpi.autodispatch.arg.broadcast.threshold = 1e8)
  inline <- npksum(y ~ x + f, data = dat, bws = c(.2, .3))
  expect_identical(inline$ksum, fit$ksum)
  expect_identical(sum(inline$timing.profile$comm_notes == "mpi.bcast.Robj2slave"), 0L)
})

test_that("npksum publication boundaries and adjacent families retain their owners", {
  skip_if_not(spawn_mpi_slaves())
  on.exit(close_mpi_slaves(), add = TRUE)
  old <- options(np.messages = FALSE, npRmpi.profile.level = "detailed",
    npRmpi.autodispatch.arg.broadcast.threshold =
      getOption("npRmpi.autodispatch.arg.broadcast.threshold", 4096L),
    npRmpi.autodispatch.arg.broadcast.threshold.regression =
      getOption("npRmpi.autodispatch.arg.broadcast.threshold.regression", 32768L))
  on.exit(options(old), add = TRUE)
  set.seed(835L)
  x <- data.frame(x = runif(80))
  y <- rnorm(80)
  size <- as.numeric(object.size(y))
  values <- list()
  for (offset in c(-1, 0, 1)) {
    options(npRmpi.autodispatch.arg.broadcast.threshold = size + offset)
    result <- npksum(txdat = x, tydat = y, bws = .2)
    values[[as.character(offset)]] <- result$ksum
    expect_identical(sum(result$timing.profile$comm_notes == "mpi.bcast.Robj2slave"), 1L)
  }
  expect_identical(values[[1]], values[[2]])
  expect_identical(values[[1]], values[[3]])
  # Non-npksum publications must not acquire the npksum bundle policy.
  options(npRmpi.autodispatch.arg.broadcast.threshold = 1L,
          npRmpi.autodispatch.arg.broadcast.threshold.regression = 1L)
  fit <- npreg(txdat = x, tydat = y, bws = .2, bandwidth.compute = FALSE)
  expect_true(all(is.finite(fitted(fit))))
  expect_true(sum(fit$timing.profile$comm_notes == "mpi.bcast.Robj2slave") > 1L)
})
