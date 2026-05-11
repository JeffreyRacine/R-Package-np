test_that("npconmode proper helper enforces binary complement probabilities", {
  helper <- getFromNamespace(".npConmodeProperProbabilities", "np")
  raw <- matrix(c(-0.2, 1.2,
                  0.4, 0.8,
                  0.3, 0.3), ncol = 2L, byrow = TRUE)

  out <- helper(raw, levels = c("0", "1"), proper = TRUE)

  expect_true(isTRUE(out$proper.requested))
  expect_true(isTRUE(out$proper.applied))
  expect_equal(rowSums(out$probabilities), rep(1, nrow(raw)), tolerance = 1e-12)
  expect_true(all(out$probabilities >= -1e-12))
  expect_equal(out$probabilities[, 2L], 1 - out$probabilities[, 1L], tolerance = 1e-12)
  expect_identical(out$proper.info$reason, "projected")
})

test_that("npconmode proper defaults follow the canonical regression type", {
  effective <- getFromNamespace(".npConmodeEffectiveProper", "np")

  expect_false(effective(list(regtype = "lc"), NULL))
  expect_true(effective(list(regtype = "ll"), NULL))
  expect_true(effective(list(regtype = "lp"), NULL))
  expect_false(effective(list(regtype = "lp", regtype.engine = "lc"), NULL))
  expect_true(effective(list(regtype = "lc", regtype.engine = "lp"), NULL))
  expect_false(effective(list(regtype = "lp"), FALSE))
  expect_true(effective(list(regtype = "lc"), TRUE))
})

test_that("npconmode formula route forwards NOMAD shortcut controls", {
  skip_if_not_installed("crs")

  set.seed(20260511)
  n <- 30L
  d <- data.frame(
    x = seq(-1, 1, length.out = n)
  )
  d$y <- factor(rbinom(n, 1L, plogis(1.5 * d$x)))

  capture.output(
    fit <- npconmode(
      y ~ x,
      data = d,
      nomad = TRUE,
      nmulti = 1L,
      nomad.nmulti = 1L,
      degree.max = 1L
    )
  )

  expect_s3_class(fit, "conmode")
  expect_true(isTRUE(fit$bws$nomad.shortcut$enabled))
  expect_identical(as.character(fit$bws$regtype.engine), "lp")
})

test_that("npconmode NOMAD fit metadata propagates through evaluation and bandwidth plots", {
  skip_if_not_installed("crs")

  count_namespace_calls <- function(fun, expr) {
    ns <- asNamespace("np")
    ctr <- new.env(parent = emptyenv())
    ctr$n <- 0L
    trace(
      what = fun,
      where = ns,
      tracer = bquote(.(ctr)$n <- .(ctr)$n + 1L),
      print = FALSE
    )
    on.exit(try(untrace(fun, where = ns), silent = TRUE), add = TRUE)
    force(expr)
    ctr$n
  }

  set.seed(20260511)
  n <- 36L
  d <- data.frame(x = seq(-1, 1, length.out = n))
  d$y <- factor(rbinom(n, 1L, plogis(1.5 * d$x)))
  nd <- data.frame(x = c(-0.75, -0.1, 0.35, 0.8))

  capture.output(
    fit <- npconmode(
      y ~ x,
      data = d,
      newdata = nd,
      nomad = TRUE,
      nmulti = 1L,
      nomad.nmulti = 1L,
      degree.max = 1L,
      probabilities = TRUE
    )
  )

  expect_s3_class(fit, "conmode")
  expect_true(isTRUE(fit$bws$nomad.shortcut$enabled))
  expect_identical(as.character(fit$bws$regtype.engine), "lp")
  expect_length(fit$conmode, nrow(nd))
  expect_length(fitted(fit), nrow(nd))
  expect_equal(nrow(fit$probabilities), nrow(nd))
  expect_true(all(fit$probabilities >= -1e-12))
  expect_equal(rowSums(fit$probabilities), rep(1, nrow(nd)), tolerance = 1e-12)

  bw.calls <- count_namespace_calls("npcdensbw", {
    out <- plot(
      fit$bws,
      xdat = d["x"],
      ydat = d["y"],
      output = "data",
      perspective = FALSE,
      neval = 6L
    )
    expect_type(out, "list")
  })
  expect_identical(bw.calls, 0L)
})
