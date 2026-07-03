test_that("genBwSelStr suppresses unknown and zero fast counts but prints known positive ones", {
  gen_bw_sel <- getFromNamespace("genBwSelStr", "np")

  s.na <- gen_bw_sel(list(num.feval = 26L, num.feval.fast = NA_integer_))
  s.zero <- gen_bw_sel(list(num.feval = 26L, num.feval.fast = 0L))
  s.pos <- gen_bw_sel(list(num.feval = 26L, num.feval.fast = 3L))

  expect_match(s.na, "Number of Function Evaluations: 26")
  expect_false(grepl("Fast CV route:", s.na, fixed = TRUE))
  expect_false(grepl("Fast CV route:", s.zero, fixed = TRUE))
  expect_true(grepl(
    "Fast CV route: 3 of 26 function evaluations",
    s.pos,
    fixed = TRUE
  ))
})

test_that("genBwSelStr prints guarded counts alongside split fast counts", {
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

  expect_true(grepl(
    "Fast CV route: 3 of 26 function evaluations",
    s.guarded,
    fixed = TRUE
  ))
  expect_true(grepl(
    "Guarded evaluations: 5",
    s.guarded,
    fixed = TRUE
  ))
  expect_true(grepl(
    "Fast CV route: 3 of 26 function evaluations",
    s.guard.na,
    fixed = TRUE
  ))
  expect_false(grepl("Guarded evaluations:", s.guard.na, fixed = TRUE))
})

test_that("genBwSelStr splits objective-cache hits from remaining fast-route savings", {
  gen_bw_sel <- getFromNamespace("genBwSelStr", "np")

  s.cache.only <- gen_bw_sel(list(
    num.feval = 26L,
    num.feval.fast = 3L,
    nn.cache = c(objective.hits = 3L, objective.visits = 10L)
  ))
  s.cache.extra <- gen_bw_sel(list(
    num.feval = 26L,
    num.feval.fast = 5L,
    nn.cache = c(objective.hits = 3L, objective.visits = 10L)
  ))

  expect_true(grepl(
    "Powell cache: 3 repeated objective lookups avoided out of 10",
    s.cache.only,
    fixed = TRUE
  ))
  expect_false(grepl("Fast CV route:", s.cache.only, fixed = TRUE))
  expect_true(grepl(
    "Powell cache: 3 repeated objective lookups avoided out of 10",
    s.cache.extra,
    fixed = TRUE
  ))
  expect_true(grepl(
    "Fast CV route: 2 of 26 function evaluations",
    s.cache.extra,
    fixed = TRUE
  ))
})
