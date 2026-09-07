index_boundary_fixture <- function(engine = "nomad", floor = 0) {
  set.seed(20260401)
  n <- 18L
  x <- sort(runif(n)); z <- sort(runif(n))
  invisible(sin(2 * pi * x) + rnorm(n, sd = .05))
  invisible(rbeta(n, 1 + x, 2 - x))
  invisible(1 + .5 * x + cos(2 * pi * z) + rnorm(n, sd = .05))
  invisible((1 + z^2) * x + rnorm(n, sd = .05))
  y <- sin(x + .5 * z) + rnorm(n, sd = .05)
  npindexbw(xdat = data.frame(x1 = x, x2 = z), ydat = y,
    bws = c(1, .5, .35), method = "ichimura", regtype = "lp",
    bernstein.basis = TRUE, degree.select = "coordinate",
    search.engine = engine, degree.min = 0L, degree.max = 1L,
    degree.verify = FALSE, bwtype = "fixed",
    scale.factor.search.lower = floor, nmulti = 2L)
}

test_that("index trial rejection does not weaken strict bandwidth validation", {
  for (h in c(0, -1, Inf, NaN)) {
    expect_null(.npindex_finalize_bandwidth(h, "fixed", 18, strict = FALSE))
    expect_error(.npindex_finalize_bandwidth(h, "fixed", 18),
                 "positive and finite")
  }
  expect_null(.npindex_finalize_bandwidth(.1, "fixed", 18, lower = .2,
                                         strict = FALSE))
  expect_error(.npindex_finalize_bandwidth(.1, "fixed", 18, lower = .2),
               "below the continuous")
  expect_identical(.npindex_finalize_bandwidth(.2, "fixed", 18, lower = .2,
                                              strict = FALSE), .2)
  expect_error(.npindex_finalize_bandwidth(0, "generalized_nn", 18), "integer")
  expect_identical(.npindex_finalize_bandwidth(4, "adaptive_nn", 18), 4)
})

test_that("index NOMAD rejects zero trials before scientific evaluation", {
  skip_if_not_installed("crs")
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  original <- .np_nomad_search
  objective <- .npindexbw_eval_objective
  rejected <- 0L
  evaluations <- 0L
  local_mocked_bindings(
    .npindexbw_eval_objective = function(...) {
      evaluations <<- evaluations + 1L
      objective(...)
    },
    .np_nomad_search = function(...) {
      args <- list(...)
      point <- args$x0
      point[length(point) - 1L] <- 0
      before <- evaluations
      result <- args$eval_fun(point)
      expect_identical(result$admissible, FALSE)
      expect_identical(result$objective, Inf)
      expect_identical(evaluations, before)
      expect_identical(result$num.feval, 0)
      rejected <<- rejected + 1L
      do.call(original, args)
    }, .package = "npRmpi")
  result <- index_boundary_fixture()
  expect_gt(rejected, 0L)
  expect_true(is.finite(result$fval))
  expect_gt(result$bw, 0)
})

test_that("index NOMAD preserves unexpected callback errors", {
  skip_if_not_installed("crs")
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  local_mocked_bindings(.npindexbw_eval_objective = function(...) {
    stop("index-boundary injected implementation error", call. = FALSE)
  }, .package = "npRmpi")
  expect_error(index_boundary_fixture(),
               "index-boundary injected implementation error")
})

