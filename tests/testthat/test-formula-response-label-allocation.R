test_that("formula constructors do not build response-sized temporary labels", {
  withr::local_options(np.messages = FALSE)
  set.seed(380283)
  d <- data.frame(x = runif(1200), z = runif(1200))
  d$y <- sin(3*d$x) + d$z
  seen <- new.env(parent = emptyenv()); seen$lengths <- numeric(); seen$metadata <- numeric()
  original.r <- rbandwidth
  original.si <- sibandwidth
  original.sc <- scbandwidth
  original.update <- updateBwNameMetadata
  # Serial construction and MPI metadata restoration can each allocate the
  # temporary label. Observe both before the final formula names overwrite it.
  testthat::local_mocked_bindings(
    updateBwNameMetadata = function(nameList, bws) {
      seen$metadata <- c(seen$metadata, sum(nchar(nameList$ynames)))
      original.update(nameList, bws)
    },
    rbandwidth = function(ynames, ...) {
      seen$lengths <- c(seen$lengths, sum(nchar(ynames)))
      original.r(ynames = ynames, ...)
    },
    sibandwidth = function(ynames, ...) {
      seen$lengths <- c(seen$lengths, sum(nchar(ynames)))
      original.si(ynames = ynames, ...)
    },
    scbandwidth = function(ynames, ...) {
      seen$lengths <- c(seen$lengths, sum(nchar(ynames)))
      original.sc(ynames = ynames, ...)
    }, .package = "np")
  for (family in c("npreg", "npindex", "npscoef")) {
    seen$lengths <- numeric(); seen$metadata <- numeric()
    ctor <- get(paste0(family, "bw"))
    fo <- switch(family, npreg = y~x, npindex = y~x+z, npscoef = y~x|z)
    bw <- if (family == "npindex") c(1, .5, .3) else .3
    counter <- new.env(); counter$n <- 0L
    b <- ctor(fo, data = {counter$n <- counter$n + 1L; d},
              bws = bw, bandwidth.compute = FALSE)
    expect_identical(counter$n, 1L)
    expect_identical(b$ynames, "y")
    expect_gt(length(seen$lengths), 0L)
    expect_lt(max(seen$lengths), 100)
    expect_gt(length(seen$metadata), 0L)
    expect_lt(max(seen$metadata), 100)
  }
})
