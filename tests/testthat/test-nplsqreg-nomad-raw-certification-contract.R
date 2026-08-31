test_that("nplsqreg NOMAD payload certification recognizes only raw objectives", {
  valid <- getFromNamespace(".nplsqreg_raw_objective_valid", "np")

  expect_true(valid(0))
  expect_true(valid(-1))
  expect_false(valid(.Machine$double.xmax))
  expect_false(valid(Inf))
  expect_false(valid(NA_real_))
  expect_false(valid(numeric()))
})

test_that("nplsqreg NOMAD payload selection preserves valid direct and Powell owners", {
  select <- getFromNamespace(".nplsqreg_nomad_select_payload", "np")
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

make_nplsqreg_nomad_invalid_direct_fixture <- function() {
  x <- data.frame(x = seq(-1, 1, length.out = 24L))
  list(
    xdat = x,
    ydat = sin(3 * x$x) + seq_len(24L) / 1000,
    scale = rep.int(1, 24L),
    bws = 1e-9,
    nomad = TRUE,
    degree.min = 1L,
    degree.max = 2L,
    bwtype = "fixed",
    ckertype = "epanechnikov",
    nmulti = 1L,
    nomad.opts = list(MAX_BB_EVAL = 1L)
  )
}

test_that("nplsqreg NOMAD-only rejects an invalid terminal payload", {
  arguments <- make_nplsqreg_nomad_invalid_direct_fixture()
  arguments$search.engine <- "nomad"
  set.seed(937L)

  expect_error(
    do.call(nplsqregbw, arguments),
    "nplsqregbw NOMAD route did not return a raw-valid solution",
    fixed = TRUE
  )
})

test_that("nplsqreg Powell can recover an invalid direct NOMAD payload", {
  arguments <- make_nplsqreg_nomad_invalid_direct_fixture()
  arguments$search.engine <- "nomad+powell"
  arguments$itmax <- 30L
  arguments$powell.remin <- FALSE
  set.seed(938L)

  bw <- do.call(nplsqregbw, arguments)

  expect_s3_class(bw, "lsqregressionbandwidth")
  expect_true(is.finite(bw$objective))
  expect_false(identical(bw$objective, .Machine$double.xmax))
  expect_equal(bw$objective, 0.01174402187, tolerance = 1e-9)
  expect_equal(bw$bw, 0.07777644704, tolerance = 1e-9)
  expect_identical(bw$degree, 1L)
  expect_gte(bw$num.feval, 114)
  expect_identical(bw$reg.bws$nomad.restart.fval, 3.374004)
})
