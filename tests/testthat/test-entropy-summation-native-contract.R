test_that("entropy fast-route policy excludes every non-default estimator", {
  route <- getFromNamespace(
    ".np_entropy_uses_default_fixed_gaussian", "npRmpi"
  )

  expect_true(route(list()))
  expect_true(route(list(bwtype = "fixed")))
  expect_true(route(list(ckertype = "gaussian", ckerorder = 2L)))
  expect_true(route(list(
    bwtype = "fixed", ckertype = "gaussian", ckerorder = 2
  )))
  expect_false(route(list(bwtype = "generalized_nn")))
  expect_false(route(list(bwtype = "adaptive_nn")))
  expect_false(route(list(ckertype = "epanechnikov")))
  expect_false(route(list(ckerorder = 4L)))
  expect_false(route(list(bwtype = "fixed", irrelevant = TRUE)))
})

test_that("entropy count chunks stay within the fixed workspace budget", {
  chunk.size <- getFromNamespace(".np_entropy_count_chunk_size", "npRmpi")

  expect_identical(chunk.size(100L, 80, 16L), 16L)
  expect_identical(chunk.size(200000L, 80, 16L), 1L)
  expect_identical(chunk.size(100L, 8, 64L), 64L)
  expect_error(
    chunk.size(0L, 80, 16L),
    "invalid entropy count-chunk dimensions"
  )
})

test_that("entropy count plans preserve tsboot draws and RNG state", {
  drawer.factory <- getFromNamespace(
    ".np_entropy_tsboot_counts_drawer", "npRmpi"
  )
  n <- 17L
  n.sim <- 11L
  B <- 23L

  for (sim in c("fixed", "geom")) {
    blocklen <- if (identical(sim, "fixed")) 1L else 4L
    set.seed(7710L + blocklen)
    drawer <- drawer.factory(
      n = n, B = B, blocklen = blocklen, sim = sim, n.sim = n.sim
    )
    counts <- cbind(drawer(1L, 7L), drawer(8L, B))
    candidate.seed <- .Random.seed

    set.seed(7710L + blocklen)
    reference <- boot::tsboot(
      tseries = seq_len(n),
      statistic = identity,
      R = B,
      n.sim = n.sim,
      l = blocklen,
      sim = sim,
      orig.t = FALSE
    )$t
    reference.seed <- .Random.seed
    reference.counts <- vapply(
      seq_len(B),
      function(replication) tabulate(
        as.integer(reference[replication, ]), nbins = n
      ),
      integer(n)
    )

    expect_equal(counts, reference.counts, tolerance = 0)
    expect_identical(candidate.seed, reference.seed)
  }
})

test_that("native symmetry counts preserve duplicate-sample statistics", {
  statistic <- getFromNamespace(
    ".np_entropy_symmetric_gaussian_summation", "npRmpi"
  )
  support <- c(-2.3, -1.1, -0.4, 0.2, 0.7, 1.8, 2.9)
  counts <- cbind(
    c(0, 2, 1, 3, 0, 1, 2),
    c(1, 0, 4, 0, 2, 1, 1)
  )
  bandwidth <- 0.61
  candidate <- .Call(
    "C_np_entropy_symmetric_summation_counts",
    support,
    matrix(as.double(counts), nrow = length(support)),
    bandwidth,
    PACKAGE = "npRmpi"
  )
  reference <- vapply(
    seq_len(ncol(counts)),
    function(column) statistic(rep(support, counts[, column]), bandwidth),
    numeric(1L)
  )

  expect_equal(candidate, reference,
               tolerance = 128 * .Machine$double.eps)
})

test_that("native univariate counts preserve duplicate-sample statistics", {
  statistic <- getFromNamespace(
    ".np_entropy_univariate_gaussian_summation", "npRmpi"
  )
  support <- c(-2, -0.8, -0.1, 0.6, 1.5, 2.4)
  counts.x <- cbind(
    c(2, 0, 3, 1, 0, 2),
    c(0, 4, 0, 2, 1, 1)
  )
  counts.y <- cbind(
    c(0, 2, 1, 0, 3, 1),
    c(3, 0, 1, 2, 0, 1)
  )
  bandwidths <- c(0.52, 0.67)
  candidate <- .Call(
    "C_np_entropy_univariate_summation_counts",
    support,
    matrix(as.double(counts.x), nrow = length(support)),
    matrix(as.double(counts.y), nrow = length(support)),
    bandwidths,
    PACKAGE = "npRmpi"
  )
  reference <- vapply(
    seq_len(ncol(counts.x)),
    function(column) statistic(
      rep(support, counts.x[, column]),
      rep(support, counts.y[, column]),
      bandwidths[1L], bandwidths[2L]
    ),
    numeric(1L)
  )

  expect_equal(candidate, reference,
               tolerance = 128 * .Machine$double.eps)
})

test_that("explicit nearest-neighbour entropy routes bypass fixed counts", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  namespace <- asNamespace("npRmpi")
  trace(
    ".np_entropy_univariate_iid_summation_bootstrap",
    tracer = quote(stop("fixed count route entered")),
    print = FALSE,
    where = namespace
  )
  trace(
    ".np_entropy_symmetric_summation_bootstrap",
    tracer = quote(stop("fixed symmetry count route entered")),
    print = FALSE,
    where = namespace
  )
  on.exit({
    untrace(
      ".np_entropy_univariate_iid_summation_bootstrap",
      where = namespace
    )
    untrace(
      ".np_entropy_symmetric_summation_bootstrap",
      where = namespace
    )
  }, add = TRUE)

  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  x <- c(-1.8, -1.1, -0.5, -0.1, 0.3, 0.8, 1.4, 2.2)
  y <- c(-1.4, -0.7, -0.2, 0.4, 1.0, 1.7)

  for (bwtype in c("generalized_nn", "adaptive_nn")) {
    expect_s3_class(npunitest(
      x, y, method = "summation", boot.num = 9L,
      bw.x = 3L, bw.y = 3L, bwtype = bwtype,
      random.seed = 9182L
    ), "unitest")
    expect_s3_class(npsymtest(
      x, method = "summation", boot.num = 9L, bw = 3L,
      bwtype = bwtype, random.seed = 2819L
    ), "symtest")
  }
})

test_that("fast symmetry bootstrap retains its historical matrix shape", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  out <- npsymtest(
    c(-2, -1.1, -0.4, 0, 0.3, 0.9, 1.7, 2.8),
    method = "summation", boot.num = 9L, bw = 0.58,
    random.seed = 6193L
  )

  expect_identical(dim(out$Srho.bootstrap), c(9L, 1L))
})
