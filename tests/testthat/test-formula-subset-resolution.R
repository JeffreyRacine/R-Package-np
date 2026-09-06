test_that("formula subsets select rows once with data-column precedence", {
  skip_on_cran()
  expect_true(spawn_mpi_slaves(1L))
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  withr::local_preserve_seed()
  withr::local_options(np.messages = FALSE)
  d <- data.frame(x = seq(-1, 1, length.out = 12L),
                  y = cos(seq_len(12L)), keep = rep(c(TRUE, FALSE), 6L))
  keep <- rep(FALSE, 12L)
  data.count <- subset.count <- 0L
  actual <- npksum(y ~ x,
    data = { data.count <- data.count + 1L; d },
    subset = { subset.count <<- subset.count + 1L; keep }, bws = .4)
  expected <- npksum(txdat = d[d$keep, "x", drop = FALSE],
                     tydat = d$y[d$keep], bws = .4)
  expect_equal(actual$ksum, expected$ksum, tolerance = 1e-12)
  expect_error(npksum(y ~ x, data = d, subset = stop("subset sentinel"),
                     bws = .4), "subset sentinel")
  expect_identical(data.count, 1L)
  expect_identical(subset.count, 1L)
  selections <- list(NULL, c(1L, 3L, 7L),
                     c(TRUE, NA, FALSE, rep(TRUE, 9L)))
  for (selected in selections) {
    rows <- if (is.null(selected)) d else stats::na.omit(d[selected, ])
    actual <- npksum(y ~ x, data = d, subset = selected, bws = .4)
    expected <- npksum(txdat = rows["x"], tydat = rows$y, bws = .4)
    expect_equal(actual$ksum, expected$ksum, tolerance = 1e-12)
  }
  actual <- npksum(y ~ x, data = d, bws = .4)
  expected <- npksum(txdat = d["x"], tydat = d$y, bws = .4)
  expect_equal(actual$ksum, expected$ksum, tolerance = 0)
  actual <- npksum(y ~ x, data = d, subset = keep, bws = 3,
                   bwtype = "adaptive_nn")
  expected <- npksum(txdat = d[d$keep, "x", drop = FALSE],
                     tydat = d$y[d$keep], bws = 3, bwtype = "adaptive_nn")
  expect_equal(actual$ksum, expected$ksum, tolerance = 1e-12)
  # The statistic must use exactly the selected native rows.
  withr::local_preserve_seed()
  withr::local_options(np.messages = FALSE)
  set.seed(1L)
  d <- data.frame(y = rnorm(60L), x = rnorm(60L))
  d$y <- d$x + .3 * rnorm(60L)
  keep <- seq_len(nrow(d)) <= 40L
  m <- lm(y ~ x, data = d[keep, ], x = TRUE, y = TRUE)
  actual <- npcmstest(y ~ x, data = d, subset = keep, model = m,
    B = 19L, nmulti = 1L, random.seed = 42L)
  expected <- npcmstest(xdat = d[keep, "x", drop = FALSE], ydat = d$y[keep],
    model = m, B = 19L, nmulti = 1L, random.seed = 42L)
  expect_equal(actual[c("Jn", "In", "P")], expected[c("Jn", "In", "P")],
    tolerance = 1e-12)
  # Discriminate the repaired data-mask contract from the old caller-only
  # resolution: the two bindings select different rows.
  d$keep <- keep
  keep <- !keep
  shadowed <- npcmstest(y ~ x, data = d, subset = keep, model = m,
    B = 19L, nmulti = 1L, random.seed = 42L)
  expect_equal(shadowed[c("Jn", "In", "P")], expected[c("Jn", "In", "P")],
    tolerance = 1e-12)
})
