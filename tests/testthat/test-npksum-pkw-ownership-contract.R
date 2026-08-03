pkw_source_file <- function() {
  candidates <- c(
    testthat::test_path("..", "..", "src", "jksum.c"),
    testthat::test_path("..", "..", "..", "src", "jksum.c"),
    file.path(Sys.getenv("R_PACKAGE_SOURCE", ""), "src", "jksum.c"),
    file.path(getwd(), "src", "jksum.c"),
    file.path(getwd(), "..", "src", "jksum.c")
  )
  candidates <- unique(candidates[nzchar(candidates)])
  hits <- candidates[file.exists(candidates)]
  if (!length(hits)) NULL else hits[[1L]]
}

test_that("derivative kernel-weight ownership is invocation scoped", {
  source <- pkw_source_file()
  skip_if(is.null(source), "source file src/jksum.c unavailable")
  text <- paste(readLines(source, warn = FALSE), collapse = "\n")

  expect_true(grepl("typedef struct {", text, fixed = TRUE))
  expect_true(grepl("} NPPermutationWeightOutput;", text, fixed = TRUE))
  expect_match(text, "NPPermutationWeightOutput \\* const pkw_output")
  expect_false(grepl("kernel_weighted_sum_pkw_extern", text, fixed = TRUE))
  expect_false(grepl("kernel_weighted_sum_pkw_nvar_extern", text, fixed = TRUE))
  expect_false(grepl("kernel_weighted_sum_pkw_sparse_extern", text, fixed = TRUE))
  expect_false(grepl("np_jksum_tree_pkw_is_sparse", text, fixed = TRUE))
})

test_that("a rejected derivative request does not poison the next p.kw call", {
  if ("npRmpi" %in% loadedNamespaces() &&
      exists("spawn_mpi_slaves", mode = "function") &&
      !spawn_mpi_slaves())
    skip("Could not spawn MPI slaves")

  old <- options(np.tree = FALSE, np.messages = FALSE)
  on.exit(options(old), add = TRUE)

  x <- data.frame(
    x1 = seq(0.05, 0.95, length.out = 31L),
    x2 = seq(0.95, 0.05, length.out = 31L)
  )
  common <- list(
    bws = c(0.21, 0.29),
    txdat = x,
    exdat = x[seq.int(1L, nrow(x), by = 2L), , drop = FALSE],
    ckertype = "epanechnikov",
    ckerorder = 4L,
    return.kernel.weights = TRUE,
    return.derivative.kernel.weights = TRUE,
    permutation.operator = "derivative"
  )

  expect_error(
    do.call(npksum, c(common, list(compute.score = TRUE))),
    "compute.score cannot be combined with a permutation operator",
    fixed = TRUE
  )

  dense <- do.call(npksum, common)
  options(np.tree = TRUE)
  tree <- do.call(npksum, common)
  repeated <- do.call(npksum, common)

  expect_identical(tree$p.kw, dense$p.kw)
  expect_identical(repeated$p.kw, tree$p.kw)
  expect_true(all(is.finite(tree$p.kw)))

  adaptive_common <- common
  adaptive_common$bws <- c(9L, 11L)
  adaptive_common$bwtype <- "adaptive_nn"
  options(np.tree = FALSE)
  adaptive_dense <- do.call(npksum, adaptive_common)
  options(np.tree = TRUE)
  adaptive_tree <- do.call(npksum, adaptive_common)
  expect_identical(adaptive_tree$p.kw, adaptive_dense$p.kw)
  expect_identical(adaptive_tree$kw, adaptive_dense$kw)
  expect_identical(adaptive_tree$p.ksum, adaptive_dense$p.ksum)
})

test_that("mixed derivative blocks use their own tree support", {
  if ("npRmpi" %in% loadedNamespaces() &&
      exists("spawn_mpi_slaves", mode = "function") &&
      !spawn_mpi_slaves())
    skip("Could not spawn MPI slaves")

  old <- options(np.tree = FALSE, np.messages = FALSE)
  on.exit(options(old), add = TRUE)

  n <- 27L
  x <- data.frame(
    x = seq(0.04, 0.96, length.out = n),
    u = factor(rep(c("a", "b", "c"), length.out = n)),
    o = ordered(
      rep(c("low", "middle", "high"), length.out = n),
      levels = c("low", "middle", "high")
    )
  )
  bws <- c(0.18, 0.23, 0.31)
  common <- list(
    txdat = x,
    exdat = x,
    bws = bws,
    ckertype = "epanechnikov",
    ckerorder = 4L,
    permutation.operator = "derivative",
    compute.ocg = TRUE,
    return.kernel.weights = TRUE,
    return.derivative.kernel.weights = TRUE
  )

  dense <- do.call(npksum, common)
  options(np.tree = TRUE)
  tree <- do.call(npksum, common)
  repeated <- do.call(npksum, common)

  options(np.tree = FALSE)
  derivative <- npksum(
    txdat = x, exdat = x, bws = bws,
    ckertype = "epanechnikov", ckerorder = 4L,
    permutation.operator = "derivative",
    return.kernel.weights = TRUE,
    return.derivative.kernel.weights = TRUE
  )
  unordered_eval <- x
  unordered_eval$u <- factor(
    rep(levels(x$u)[1L], n), levels = levels(x$u)
  )
  unordered_reference <- npksum(
    txdat = x, exdat = unordered_eval, bws = bws,
    ckertype = "epanechnikov", ckerorder = 4L,
    return.kernel.weights = TRUE
  )
  ordered_eval <- x
  ordered_index <- ifelse(
    as.integer(x$o) == 1L, 2L, as.integer(x$o) - 1L
  )
  ordered_eval$o <- ordered(
    levels(x$o)[ordered_index], levels = levels(x$o)
  )
  ordered_reference <- npksum(
    txdat = x, exdat = ordered_eval, bws = bws,
    ckertype = "epanechnikov", ckerorder = 4L,
    return.kernel.weights = TRUE
  )

  expect_identical(tree$p.kw, dense$p.kw)
  expect_identical(repeated$p.kw, tree$p.kw)
  expect_identical(dense$p.kw[, , 1L], derivative$p.kw)
  expect_identical(dense$p.kw[, , 2L], unordered_reference$kw)
  expect_identical(dense$p.kw[, , 3L], ordered_reference$kw)
  expect_true(all(is.finite(tree$p.kw)))
})
