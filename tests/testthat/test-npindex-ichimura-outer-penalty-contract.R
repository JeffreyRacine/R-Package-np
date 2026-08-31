test_that("Ichimura outer guidance uses the CVLS null contract", {
  penalty <- getFromNamespace(".npindexbw_ichimura_outer_penalty", "npRmpi")
  y <- c(-3, -1, 2, 6)
  n <- length(y)
  reference <- (n / (n - 1))^2 * mean((y - mean(y))^2)

  expect_identical(penalty(y), 10 * reference)
  expect_identical(penalty(y + 1024), penalty(y))
  expect_identical(penalty(-2 * y), 4 * penalty(y))

  constant.penalty <- penalty(rep(7, 8))
  expect_identical(
    constant.penalty,
    .Machine$double.xmin * .Machine$double.eps
  )
  expect_gt(constant.penalty, 0)
  expect_identical(penalty(1), .Machine$double.xmax)
  expect_identical(penalty(c(1, NA_real_)), .Machine$double.xmax)
  expect_identical(penalty(c(-1e308, 1e308)), .Machine$double.xmax)
})

test_that("Ichimura outer mapping changes only exact raw invalidity", {
  map <- getFromNamespace(".npindexbw_map_ichimura_outer_result", "npRmpi")
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

test_that("Ichimura regression leaf requests raw-or-DBL_MAX and propagates errors", {
  leaf <- getFromNamespace(".npindexbw_eval_ichimura_lp_via_npreg", "npRmpi")
  mode.seen <- NA_character_
  build.leaf <- function(...) {
    list(
      xdat = data.frame(index = 1:4),
      bws = structure(list(), class = "rbandwidth")
    )
  }
  args <- list(
    index = 1:4,
    ydat = c(1, 2, 3, 4),
    h = 1,
    bws = structure(list(), class = "sibandwidth"),
    spec = list()
  )

  out <- testthat::with_mocked_bindings(
    .npindexbw_lp_regression_leaf = build.leaf,
    .np_progress_with_nested_bandwidth_heartbeat = function(expr) force(expr),
    .npregbw_eval_only = function(..., invalid.penalty, penalty.multiplier) {
      mode.seen <<- invalid.penalty
      list(objective = .Machine$double.xmax, num.feval.fast = 0)
    },
    .package = "npRmpi",
    do.call(leaf, args)
  )
  expect_identical(mode.seen, "dbmax")
  expect_identical(out$objective, .Machine$double.xmax)

  expect_error(
    testthat::with_mocked_bindings(
      .npindexbw_lp_regression_leaf = build.leaf,
      .np_progress_with_nested_bandwidth_heartbeat = function(expr) force(expr),
      .npregbw_eval_only = function(...) stop("implementation defect", call. = FALSE),
      .package = "npRmpi",
      do.call(leaf, args)
    ),
    "implementation defect",
    fixed = TRUE
  )
})
