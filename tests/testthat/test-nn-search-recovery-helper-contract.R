test_that("ordinary NN recovery schedules only declared coordinates", {
  point <- c(2, 0.375, 3, 0.625, 7)
  schedule <- np:::.np_nn_ordinary_schedule(
    point = point,
    nn.indices = c(1L, 3L),
    caps = c(9L, 6L)
  )

  expect_equal(
    lapply(schedule, function(x) x[c(1L, 3L)]),
    list(c(4, 6), c(8, 6), c(9, 6))
  )
  for (candidate in schedule)
    expect_identical(candidate[c(2L, 4L, 5L)], point[c(2L, 4L, 5L)])
  expect_identical(schedule[[length(schedule)]][c(1L, 3L)], c(9, 6))
})

test_that("ordinary NN recovery uses raw validity and one deterministic order", {
  visited <- list()
  raw.eval <- function(point) {
    visited[[length(visited) + 1L]] <<- point
    if (point[[1L]] < 8) {
      np:::.np_nn_abort_candidate_invalid(
        "synthetic raw candidate is invalid",
        owner = "synthetic",
        point = point
      )
    }
    if (point[[3L]] < 7) return(.Machine$double.xmax)
    1.25
  }

  result <- np:::.np_nn_find_raw_valid_start(
    point = c(2, 0.3, 2, 0.7),
    nn.indices = c(1L, 3L),
    caps = c(9L, 7L),
    raw.eval = raw.eval
  )

  expect_true(result$found)
  expect_identical(result$point, c(8, 0.3, 7, 0.7))
  expect_identical(result$objective, 1.25)
  expect_identical(result$evaluations, 2L)
  expect_identical(visited, list(c(4, 0.3, 4, 0.7), c(8, 0.3, 7, 0.7)))
})

test_that("ordinary NN recovery propagates unknown errors", {
  expect_error(
    np:::.np_nn_find_raw_valid_start(
      point = c(2, 2),
      nn.indices = 1:2,
      caps = c(8, 8),
      raw.eval = function(point) stop("unknown synthetic defect")
    ),
    "unknown synthetic defect",
    fixed = TRUE
  )
})

test_that("extended NN incumbents skip ordinary recovery without clamping", {
  evaluated <- FALSE
  result <- np:::.np_nn_find_raw_valid_start(
    point = c(12, 2),
    nn.indices = 1:2,
    caps = c(9, 9),
    raw.eval = function(point) {
      evaluated <<- TRUE
      1
    }
  )

  expect_false(result$found)
  expect_null(result$point)
  expect_identical(result$evaluations, 0L)
  expect_false(evaluated)
  expect_error(
    np:::.np_nn_find_raw_valid_start(
      point = c(0, 2), nn.indices = 1:2, caps = c(9, 9),
      raw.eval = function(point) 1
    ),
    "outside the ordinary NN domain",
    fixed = TRUE
  )
})

test_that("raw objective certification rejects penalties at the terminal scale", {
  expect_true(np:::.np_nn_raw_objective_valid(0))
  expect_true(np:::.np_nn_raw_objective_valid(-1e7))
  expect_false(np:::.np_nn_raw_objective_valid(NA_real_))
  expect_false(np:::.np_nn_raw_objective_valid(Inf))
  expect_false(np:::.np_nn_raw_objective_valid(.Machine$double.xmax))
  expect_false(np:::.np_nn_raw_objective_valid(-.Machine$double.xmax))
})

test_that("raw endpoint certification returns a scalar or a typed condition", {
  point <- c(4, 0.25)
  expect_identical(
    np:::.np_nn_certify_raw_point(
      point, raw.eval = function(x) 1.5, owner = "synthetic MADS"
    ),
    1.5
  )

  condition <- tryCatch(
    np:::.np_nn_certify_raw_point(
      point,
      raw.eval = function(x) .Machine$double.xmax,
      owner = "synthetic MADS"
    ),
    np_nn_candidate_invalid = identity
  )
  expect_s3_class(condition, "np_nn_candidate_invalid")
  expect_identical(condition$owner, "synthetic MADS")
  expect_identical(condition$point, point)
  expect_identical(condition$raw.objective, .Machine$double.xmax)
  expect_identical(
    np:::.np_nn_certify_raw_value(1.5, point, "synthetic MADS"),
    1.5
  )
  direct.condition <- tryCatch(
    np:::.np_nn_certify_raw_value(
      .Machine$double.xmax, point, "synthetic MADS"
    ),
    np_nn_candidate_invalid = identity
  )
  expect_identical(direct.condition, condition)
})
