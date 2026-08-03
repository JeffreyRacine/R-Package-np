test_that("computed generalized-NN convolution owns training bandwidths", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(), add = TRUE)

  z <- c(
    -1.80, -1.10, -0.20, 0.35, 1.05, 1.80,
    seq(-1.67, 1.67, length.out = 48L)
  )
  x <- data.frame(z = z)
  neighbour <- 7L
  bandwidth <- vapply(
    z,
    function(evaluation) sort(abs(z - evaluation), method = "radix")[[neighbour]],
    numeric(1)
  )

  native <- npksum(
    bws = neighbour,
    txdat = x,
    exdat = x,
    bwtype = "generalized_nn",
    ckertype = "gaussian",
    operator = "convolution",
    return.kernel.weights = TRUE
  )

  oracle <- outer(
    seq_along(z),
    seq_along(z),
    Vectorize(function(training, evaluation) {
      h.training <- bandwidth[[training]]
      h.evaluation <- bandwidth[[evaluation]]
      h.combined <- sqrt(h.training^2 + h.evaluation^2)

      0.3989422803 * h.training * h.evaluation / h.combined *
        exp(-0.5 * (z[[evaluation]] - z[[training]])^2 / h.combined^2)
    })
  )

  expect_equal(native$kw, oracle, tolerance = 2e-14)
  expect_equal(as.double(native$ksum), colSums(oracle), tolerance = 2e-14)
})

test_that("alternate convolution bandwidth storage has one explicit owner", {
  roots <- unique(c(
    Sys.getenv("NP_SOURCE_ROOT", unset = ""),
    normalizePath(file.path(getwd(), "..", ".."), mustWork = FALSE),
    normalizePath(getwd(), mustWork = FALSE)
  ))
  source.path <- NULL
  for (root in roots[nzchar(roots)]) {
    candidate <- file.path(root, "src", "jksum.c")
    if (file.exists(candidate)) {
      source.path <- candidate
      break
    }
  }
  skip_if(is.null(source.path), "package source unavailable")
  engine <- paste(readLines(source.path, warn = FALSE), collapse = "\n")

  expect_match(
    engine,
    "np_gnn_convolution_training_bandwidth_prepare(",
    fixed = TRUE
  )
  expect_match(
    engine,
    "kernel, BW_GEN_NN, num_obs_train, num_obs_train,",
    fixed = TRUE
  )
  expect_match(
    engine,
    paste0(
      "if(!bandwidth_provided && (BANDWIDTH_reg == BW_GEN_NN) &&\n",
      "     any_convolution && (num_reg_continuous > 0))\n",
      "    free_tmat(matrix_alt_bandwidth);"
    ),
    fixed = TRUE
  )
  expect_false(grepl(
    "BW_GEN_NN) && any_convolution))\\n    free_tmat(matrix_alt_bandwidth)",
    engine,
    fixed = TRUE
  ))
})
