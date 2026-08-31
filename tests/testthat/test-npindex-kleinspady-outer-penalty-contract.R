test_that("Klein-Spady outer guidance uses the clipped CVKS null", {
  penalty <- getFromNamespace(".npindexbw_kleinspady_outer_penalty", "np")
  reference <- function(y) {
    n <- length(y)
    ones <- sum(y == 1)
    zeros <- n - ones
    floor <- sqrt(.Machine$double.eps)
    clip <- function(p) min(1 - floor, max(floor, p))
    loss <- 0
    if (zeros > 0)
      loss <- zeros * -log1p(-clip(ones / (n - 1)))
    if (ones > 0)
      loss <- loss + ones * -log(clip((ones - 1) / (n - 1)))
    loss / n
  }

  for (y in list(
    c(0, 1), c(0, 0, 0, 1, 1), c(rep(0, 9), 1),
    c(0, rep(1, 9)), rep(0, 8), rep(1, 8)
  )) {
    expect_equal(penalty(y), 10 * reference(y), tolerance = 8 * .Machine$double.eps)
    expect_identical(penalty(rev(y)), penalty(y))
    expect_equal(penalty(1 - y), penalty(y), tolerance = 8 * .Machine$double.eps)
    expect_gt(penalty(y), reference(y))
  }

  expect_identical(penalty(1), .Machine$double.xmax)
  expect_identical(penalty(c(0, NA_real_)), .Machine$double.xmax)
  expect_identical(penalty(c(0, 0.5, 1)), .Machine$double.xmax)
})

test_that("Klein-Spady outer mapping changes only exact raw invalidity", {
  map <- getFromNamespace(".npindexbw_map_kleinspady_outer_result", "np")
  penalty <- 17

  for (value in c(-1, 0, 17, 100)) {
    result <- list(objective = value, num.feval.fast = 3)
    expect_identical(map(result, penalty), result)
  }

  terminal <- list(objective = .Machine$double.xmax, num.feval.fast = 0)
  expect_identical(map(terminal, penalty)$objective, penalty)
  expect_error(
    map(list(objective = Inf, num.feval.fast = 0), penalty),
    "non-terminal invalid objective",
    fixed = TRUE
  )
  expect_error(
    map(list(objective = NA_real_, num.feval.fast = 0), penalty),
    "non-terminal invalid objective",
    fixed = TRUE
  )
})

test_that("Klein-Spady regression leaf requests DBL_MAX and propagates errors", {
  leaf <- getFromNamespace(".npindexbw_eval_kleinspady_lp_via_npreg", "np")
  mode.seen <- NA_character_
  build.leaf <- function(...) {
    list(
      xdat = data.frame(index = 1:4),
      bws = structure(list(), class = "rbandwidth")
    )
  }
  args <- list(
    index = 1:4, ydat = c(0, 1, 0, 1), h = 1,
    bws = structure(list(), class = "sibandwidth"), spec = list()
  )

  out <- testthat::with_mocked_bindings(
    .npindexbw_lp_regression_leaf = build.leaf,
    .np_progress_with_nested_bandwidth_heartbeat = function(expr) force(expr),
    .npregbw_eval_only = function(..., invalid.penalty, penalty.multiplier) {
      mode.seen <<- invalid.penalty
      list(objective = .Machine$double.xmax, num.feval.fast = 0)
    },
    .package = "np",
    do.call(leaf, args)
  )
  expect_identical(mode.seen, "dbmax")
  expect_identical(out$objective, .Machine$double.xmax)

  expect_error(
    testthat::with_mocked_bindings(
      .npindexbw_lp_regression_leaf = build.leaf,
      .np_progress_with_nested_bandwidth_heartbeat = function(expr) force(expr),
      .npregbw_eval_only = function(...) stop("implementation defect", call. = FALSE),
      .package = "np",
      do.call(leaf, args)
    ),
    "implementation defect",
    fixed = TRUE
  )
})
