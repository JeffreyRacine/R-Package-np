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

  s.cache.only <- gen_bw_sel(structure(list(
    num.feval = 26L,
    num.feval.fast = 3L,
    nn.cache = c(objective.hits = 3L, objective.visits = 10L)
  ), class = "bandwidth"))
  s.cache.extra <- gen_bw_sel(structure(list(
    num.feval = 26L,
    num.feval.fast = 5L,
    nn.cache = c(objective.hits = 3L, objective.visits = 10L)
  ), class = "bandwidth"))

  expect_true(grepl(
    "Evaluation cache (Powell): 3 hits / 10 lookups (30.0%)",
    s.cache.only,
    fixed = TRUE
  ))
  expect_false(grepl("Fast CV route:", s.cache.only, fixed = TRUE))
  expect_true(grepl(
    "Evaluation cache (Powell): 3 hits / 10 lookups (30.0%)",
    s.cache.extra,
    fixed = TRUE
  ))
  expect_true(grepl(
    "Fast CV route: 2 of 26 function evaluations",
    s.cache.extra,
    fixed = TRUE
  ))
})

test_that("genBwSelStr attributes evaluation caches to the proven search stage", {
  gen_bw_sel <- getFromNamespace("genBwSelStr", "np")
  cache <- c(hits = 3L, visits = 10L)

  s.optim <- gen_bw_sel(structure(list(
    num.feval = 26L,
    nn.cache = cache,
    pomethod = "BFGS"
  ), class = "scbandwidth"))
  s.legacy <- gen_bw_sel(list(num.feval = 26L, nn.cache = cache))
  s.hybrid <- gen_bw_sel(structure(list(
    num.feval = 26L,
    nn.cache = cache,
    pomethod = "Nelder-Mead",
    degree.search = list(
      engine = "nomad+powell",
      restart.results = list(list(native = list(
        cache_hits = 2L,
        total_evaluations = 8L
      )))
    )
  ), class = "scbandwidth"))

  expect_true(grepl(
    "Evaluation cache (R optim: BFGS): 3 hits / 10 lookups (30.0%)",
    s.optim,
    fixed = TRUE
  ))
  expect_true(grepl(
    "Evaluation cache: 3 hits / 10 lookups (30.0%)",
    s.legacy,
    fixed = TRUE
  ))
  expect_true(grepl(
    "Evaluation cache (NOMAD): 2 hits / 8 lookups (25.0%)",
    s.hybrid,
    fixed = TRUE
  ))
  expect_true(grepl(
    "Evaluation cache (R optim: Nelder-Mead): 3 hits / 10 lookups (30.0%)",
    s.hybrid,
    fixed = TRUE
  ))
})

test_that("genTimingStr attributes nominal Powell timing to family-native R optim", {
  gen_timing <- getFromNamespace("genTimingStr", "np")

  core <- gen_timing(structure(list(
    total.time = 1,
    nomad.time = 0.6,
    powell.time = 0.4
  ), class = "bandwidth"))
  semiparametric <- gen_timing(structure(list(
    total.time = 1,
    nomad.time = 0.6,
    powell.time = 0.4,
    pomethod = "CG"
  ), class = "sibandwidth"))

  expect_true(grepl("NOMAD 0.6s, Powell 0.4s", core, fixed = TRUE))
  expect_true(grepl("NOMAD 0.6s, R optim: CG 0.4s", semiparametric, fixed = TRUE))
  expect_false(grepl("Powell 0.4s", semiparametric, fixed = TRUE))
})
