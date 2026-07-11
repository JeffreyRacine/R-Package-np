test_that("bounded kernels retain the finite-support uniform large-h limit", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  library(npRmpi)
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  set.seed(42)
  y <- runif(250)
  a <- min(y)
  b <- max(y)
  span <- b-a
  eval <- data.frame(y = c(a, a + span*1e-14, (a+b)/2,
                           b - span*1e-14, b))
  specs <- rbind(
    data.frame(kernel = "gaussian", order = c(2L, 4L, 6L, 8L)),
    data.frame(kernel = "epanechnikov", order = c(2L, 4L, 6L, 8L)),
    data.frame(kernel = "uniform", order = 2L),
    data.frame(kernel = "truncated gaussian", order = 2L)
  )

  for (i in seq_len(nrow(specs))) {
    kernel <- as.character(specs$kernel[i])
    args <- list(
      dat = data.frame(y = y),
      bws = span/1e-16,
      bandwidth.compute = FALSE,
      ckertype = kernel,
      ckerbound = "range"
    )
    if (kernel != "uniform")
      args$ckerorder <- specs$order[i]
    bw <- if (kernel == "uniform") {
      suppressWarnings(do.call(npudensbw, args))
    } else {
      do.call(npudensbw, args)
    }
    got <- npudens(tdat = data.frame(y = y), edat = eval, bws = bw)$dens
    expect_equal(
      as.numeric(got),
      rep.int(1/span, nrow(eval)),
      tolerance = 2e-12,
      info = paste(kernel, specs$order[i])
    )
  }
})

test_that("bounded Gaussian normalization has no historical width seam", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  library(npRmpi)
  old <- options(np.messages = FALSE, np.tree = FALSE, np.largeh = FALSE)
  on.exit(options(old), add = TRUE)

  set.seed(42)
  y <- runif(250)
  a <- min(y)
  b <- max(y)
  span <- b-a
  eval <- data.frame(y = c(a, a + span*1e-14, (a+b)/2, b))
  values <- lapply(c(0.999999e-5, 1e-5, 1.000001e-5), function(width) {
    bw <- npudensbw(
      dat = data.frame(y = y), bws = span/width,
      bandwidth.compute = FALSE, ckertype = "gaussian", ckerorder = 2L,
      ckerbound = "range"
    )
    npudens(tdat = data.frame(y = y), edat = eval, bws = bw)$dens
  })
  values <- do.call(cbind, values)

  expect_lt(max(values)-min(values), 5e-10)
  expect_lt(max(abs(values - 1/span)), 5e-10)
})
