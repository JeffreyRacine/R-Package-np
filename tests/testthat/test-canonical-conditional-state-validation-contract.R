canonical_conditional_state <- function(regtype.engine = "lc",
                                        degree.engine = 0L,
                                        ncon = 1L,
                                        basis.engine = "glp",
                                        bernstein.basis.engine = FALSE) {
  list(
    regtype.engine = regtype.engine,
    basis.engine = basis.engine,
    degree.engine = degree.engine,
    bernstein.basis.engine = bernstein.basis.engine,
    xncon = ncon
  )
}

test_that("canonical conditional engine states normalize predictably", {
  validate <- getFromNamespace(
    "npValidateConditionalRegEngineState",
    "np"
  )

  lc <- validate("lc", "glp", c(0, 0), FALSE, 2L, "test")
  expect_identical(
    lc,
    list(
      reg.engine = "lc",
      basis.engine = "glp",
      degree.engine = c(0L, 0L),
      bernstein.engine = FALSE,
      ncon = 2L
    )
  )

  lp0 <- validate("lp", "tensor", 0, FALSE, 1L, "test")
  expect_identical(lp0$reg.engine, "lp")
  expect_identical(lp0$degree.engine, 0L)

  lp <- validate("lp", "additive", c(1, 3), TRUE, 2L, "test")
  expect_identical(lp$degree.engine, c(1L, 3L))
  expect_true(lp$bernstein.engine)

  all.cat <- validate("lc", "glp", integer(), FALSE, 0L, "test")
  expect_identical(all.cat$degree.engine, integer())
  expect_identical(all.cat$ncon, 0L)
})

test_that("incoherent conditional engine states fail rather than normalize", {
  validate <- getFromNamespace(
    "npValidateConditionalRegEngineState",
    "np"
  )

  expect_error(validate("ll", "glp", 1L, FALSE, 1L), "exactly 'lc' or 'lp'")
  expect_error(validate("lp", "unknown", 1L, FALSE, 1L), "basis.engine")
  expect_error(validate("lp", "glp", integer(), FALSE, 1L), "degree.engine")
  expect_error(validate("lp", "glp", NA_integer_, FALSE, 1L), "degree.engine")
  expect_error(validate("lp", "glp", 1.5, FALSE, 1L), "degree.engine")
  expect_error(validate("lp", "glp", -1L, FALSE, 1L), "degree.engine")
  expect_error(validate("lp", "glp", 101L, FALSE, 1L), "degree.engine")
  expect_error(validate("lp", "glp", 1L, NA, 1L), "bernstein.basis.engine")
  expect_error(validate("lc", "glp", 1L, FALSE, 1L), "all zero")
  expect_error(validate("lc", "glp", 0L, TRUE, 1L), "must be FALSE")
  expect_error(validate("lp", "glp", integer(), FALSE, 0L),
               "requires at least one continuous predictor")
})

test_that("conditional state extraction is exact and rejects missing duplicates", {
  spec <- getFromNamespace("npConditionalRegEngineSpec", "np")

  state <- canonical_conditional_state(
    regtype.engine = "lp",
    degree.engine = 2L
  )
  expect_identical(spec(state)$degree.engine, 2L)

  partial <- state
  names(partial)[names(partial) == "regtype.engine"] <- "regtype.engine.shadow"
  expect_error(spec(partial), "missing required field 'regtype.engine'")

  duplicate <- c(state, list(regtype.engine = "lc"))
  expect_error(spec(duplicate), "duplicate field 'regtype.engine'")

  bad.lc <- state
  bad.lc$regtype.engine <- "lc"
  expect_error(spec(bad.lc), "degree.engine.*all zero")
})

test_that("canonical LP0 detection is total only over valid engine state", {
  is.lp0 <- getFromNamespace("npIsCanonicalLp0", "np")

  expect_true(is.lp0("lc", c(0L, 0L), 2L))
  expect_true(is.lp0("lp", c(0L, 0L), 2L))
  expect_false(is.lp0("lp", c(0L, 1L), 2L))
  expect_true(is.lp0("lc", integer(), 0L))
  expect_error(is.lp0("lc", 1L, 1L), "all zero")
  expect_error(is.lp0("ll", 1L, 1L), "exactly 'lc' or 'lp'")
})

test_that("mutated public conditional bandwidth state fails at consumers", {
  set.seed(20260730)
  x <- data.frame(x = stats::runif(24L))
  y <- data.frame(y = x$x + stats::rnorm(24L, sd = 0.1))
  bw <- npcdensbw(
    xdat = x,
    ydat = y,
    bws = c(0.3, 0.3),
    bandwidth.compute = FALSE
  )

  mutated <- bw
  names(mutated)[names(mutated) == "regtype.engine"] <-
    "regtype.engine.shadow"

  expect_error(
    npcdens(
      bws = mutated,
      txdat = x,
      tydat = y,
      exdat = x[1:3, , drop = FALSE],
      eydat = y[1:3, , drop = FALSE]
    ),
    "missing required field 'regtype.engine'"
  )
})
