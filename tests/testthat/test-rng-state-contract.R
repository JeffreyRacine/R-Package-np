test_that("categorical padding helper is deterministic and RNG-pristine", {
  num_not_in <- getFromNamespace("numNotIn", "np")

  set.seed(8675309)
  before <- .Random.seed

  expect_identical(num_not_in(numeric()), 0)
  expect_identical(num_not_in(c(1, 2, 3)), 0)
  expect_false(num_not_in(c(0, 1, 2)) %in% c(0, 1, 2))
  expect_false(num_not_in(c(-3, -2, -1, 0, 1, 2, 3)) %in% c(-3, -2, -1, 0, 1, 2, 3))
  expect_true(is.finite(num_not_in(c(-.Machine$double.xmax, 0, .Machine$double.xmax))))

  expect_identical(.Random.seed, before)
})

test_that("npseed rejects malformed seeds without touching R RNG", {
  npseed <- getFromNamespace("npseed", "np")

  set.seed(314159)
  before <- .Random.seed

  expect_error(npseed(NA_real_), "single finite numeric")
  expect_error(npseed(1.5), "non-negative integer")
  expect_error(npseed(.Machine$integer.max + 1), "non-negative integer")
  expect_silent(npseed(0))

  expect_identical(.Random.seed, before)
})

test_that("npseed makes C bandwidth multistart search reproducible and R RNG-pristine", {
  npseed <- getFromNamespace("npseed", "np")
  npregbw <- getFromNamespace("npregbw", "np")
  x <- seq(-1, 1, length.out = 24)
  y <- x^2 + 0.1 * x

  old_opts <- options(np.messages = FALSE)
  on.exit(options(old_opts), add = TRUE)

  run_search <- function() {
    npseed(2468)
    npregbw(
      xdat = data.frame(x = x),
      ydat = y,
      nmulti = 2L,
      regtype = "lc",
      bwmethod = "cv.ls",
      ckertype = "gaussian"
    )
  }

  set.seed(999)
  before <- .Random.seed

  b1 <- run_search()
  mid <- .Random.seed
  b2 <- run_search()
  after <- .Random.seed

  expect_identical(mid, before)
  expect_identical(after, before)
  expect_equal(as.numeric(b1$bw), as.numeric(b2$bw), tolerance = 0)
  expect_equal(as.numeric(b1$fval), as.numeric(b2$fval), tolerance = 0)
  expect_identical(as.integer(b1$num.feval), as.integer(b2$num.feval))

  if (!is.null(b1$fval.history) || !is.null(b2$fval.history)) {
    expect_equal(as.numeric(b1$fval.history), as.numeric(b2$fval.history), tolerance = 0)
  }
})
