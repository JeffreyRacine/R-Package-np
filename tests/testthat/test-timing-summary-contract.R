library(npRmpi)

quiet_eval <- function(expr) {
  value <- NULL
  capture.output(value <- force(expr))
  value
}

test_that("NOMAD plus Powell timing is carried through bandwidth and fit summaries", {
  skip_if_not_installed("crs")
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(x = sort(runif(24)))
  dat$y <- sin(2 * pi * dat$x) + rnorm(nrow(dat), sd = 0.05)

  bw <- quiet_eval(npregbw(
    y ~ x,
    data = dat,
    regtype = "lp",
    degree.select = "coordinate",
    degree.min = 0L,
    degree.max = 2L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  ))
  fit <- quiet_eval(npreg(bws = bw, txdat = data.frame(x = dat$x), tydat = dat$y))

  expect_identical(bw$degree.search$mode, "nomad+powell")
  expect_true(is.finite(bw$nomad.time))
  expect_true(is.finite(bw$powell.time))
  expect_true(is.finite(bw$total.time))
  expect_equal(as.double(bw$total.time),
               as.double(bw$nomad.time + bw$powell.time),
               tolerance = 1e-8)
  expect_equal(as.double(bw$degree.search$optim.time),
               as.double(bw$nomad.time + bw$powell.time),
               tolerance = 1e-8)
  expect_equal(fit$optim.time, bw$total.time, tolerance = 1e-8)
  expect_equal(fit$nomad.time, bw$nomad.time, tolerance = 1e-8)
  expect_equal(fit$powell.time, bw$powell.time, tolerance = 1e-8)

  bw_txt <- paste(capture.output(summary(bw)), collapse = "\n")
  fit_txt <- paste(capture.output(summary(fit)), collapse = "\n")
  expect_match(bw_txt, "NOMAD ", fixed = TRUE)
  expect_match(bw_txt, "Powell ", fixed = TRUE)
  expect_match(bw_txt, "MPI Session:", fixed = TRUE)
  expect_match(bw_txt, "MPI Call Profile:", fixed = TRUE)
  expect_match(fit_txt, "NOMAD ", fixed = TRUE)
  expect_match(fit_txt, "Powell ", fixed = TRUE)
  expect_match(fit_txt, "fit ", fixed = TRUE)
  expect_match(fit_txt, "MPI Session:", fixed = TRUE)
  expect_match(fit_txt, "MPI Call Profile:", fixed = TRUE)
})
