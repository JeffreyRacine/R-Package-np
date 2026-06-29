test_that("genBwSelStr suppresses unknown fast counts but prints known ones", {
  gen_bw_sel <- getFromNamespace("genBwSelStr", "np")

  s.na <- gen_bw_sel(list(num.feval = 26L, num.feval.fast = NA_integer_))
  s.zero <- gen_bw_sel(list(num.feval = 26L, num.feval.fast = 0L))
  s.pos <- gen_bw_sel(list(num.feval = 26L, num.feval.fast = 3L))

  expect_match(s.na, "Number of Function Evaluations: 26")
  expect_false(grepl("fast =", s.na, fixed = TRUE))
  expect_false(grepl("Fast CV route:", s.zero, fixed = TRUE))
  expect_match(s.pos, "Fast CV route:\\s+3 of 26 function evaluations")
})

test_that("genBwSelStr prints guarded counts alongside fast counts", {
  gen_bw_sel <- getFromNamespace("genBwSelStr", "np")

  s.guarded <- gen_bw_sel(list(
    num.feval = 26L,
    num.feval.fast = 3L,
    num.feval.guarded = 5L
  ))
  s.guard.na <- gen_bw_sel(list(
    num.feval = 26L,
    num.feval.fast = 3L,
    num.feval.guarded = NA_real_
  ))

  expect_match(s.guarded, "Fast CV route:\\s+3 of 26 function evaluations")
  expect_true(grepl("Guarded evaluations: 5", s.guarded, fixed = TRUE))
  expect_match(s.guard.na, "Fast CV route:\\s+3 of 26 function evaluations")
  expect_false(grepl("Guarded evaluations:", s.guard.na, fixed = TRUE))
})
