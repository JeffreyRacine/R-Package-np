test_that("npRmpi nplsqreg NOMAD certification recognizes only raw objectives", {
  valid <- getFromNamespace(".nplsqreg_raw_objective_valid", "npRmpi")

  expect_true(valid(0))
  expect_true(valid(-1))
  expect_false(valid(.Machine$double.xmax))
  expect_false(valid(Inf))
  expect_false(valid(NA_real_))
  expect_false(valid(numeric()))
})

test_that("npRmpi nplsqreg NOMAD selection preserves direct and Powell owners", {
  select <- getFromNamespace(".nplsqreg_nomad_select_payload", "npRmpi")
  direct <- list(objective = 2, owner = "direct")
  hot.better <- list(objective = 1, owner = "powell")
  hot.worse <- list(objective = 3, owner = "powell")

  expect_identical(select(direct)$payload$owner, "direct")
  expect_identical(select(direct, hot.better)$payload$owner, "powell")
  expect_identical(select(direct, hot.worse)$payload$owner, "direct")
  expect_identical(
    select(list(objective = .Machine$double.xmax), hot.better)$payload$owner,
    "powell"
  )
  expect_error(
    select(
      list(objective = .Machine$double.xmax),
      list(objective = .Machine$double.xmax)
    ),
    "nplsqregbw NOMAD route did not return a raw-valid solution",
    fixed = TRUE
  )
})
