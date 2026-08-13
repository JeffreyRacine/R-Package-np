make_nplsqreg_auto_fixture <- function(n = 20L, p = 1L, seed = 20260813L) {
  set.seed(seed)
  x <- seq(0.05, 0.95, length.out = n)
  out <- data.frame(x1 = x)
  if (p >= 2L)
    out$x2 <- cos(seq(0.1, 2.1, length.out = n))
  if (p >= 3L)
    out$x3 <- sin(seq(0.2, 2.2, length.out = n))
  out$y <- sin(2 * pi * x) + stats::rnorm(n, sd = 0.04)
  out
}

test_that("nplsqreg AUTO uses the shared one-dimensional exhaustive policy", {
  options(np.messages = FALSE, np.tree = FALSE)
  dat <- make_nplsqreg_auto_fixture()
  bw <- nplsqregbw(
    y ~ x1,
    data = dat,
    scale = rep.int(1, nrow(dat)),
    nomad = "auto",
    degree.min = 0L,
    degree.max = 1L,
    nmulti = 1L,
    optim.control = list(maxit = 2L)
  )

  expect_s3_class(bw, "lsqregressionbandwidth")
  expect_identical(bw$reg.bws$degree.search$mode, "exhaustive")
  expect_true(isTRUE(bw$reg.bws$degree.search$certified))
  expect_identical(bw$reg.bws$nomad.shortcut$source, "auto")
  expect_identical(bw$reg.bws$nomad.shortcut$nomad, "auto")
  expect_true(all(c(0L, 1L) %in% bw$reg.bws$degree.search$trace$degree))
  expect_true(is.finite(bw$delta))
  expect_true(is.finite(bw$objective))
})
test_that("nplsqreg AUTO keeps the joint fixed-degree core in every cell", {
  source <- paste(deparse(getFromNamespace(".nplsqreg_cell_search", "np")),
                  collapse = "\n")
  expect_match(source, ".nplsqreg_call_fixed_degree_core", fixed = TRUE)
  expect_match(source, "bandwidth.compute = TRUE", fixed = TRUE)
  expect_match(source, ".np_degree_search", fixed = TRUE)
  expect_match(source, "direction = \"min\"", fixed = TRUE)
})

test_that("nplsqreg AUTO preserves every supported bandwidth type", {
  options(np.messages = FALSE, np.tree = FALSE)
  dat <- make_nplsqreg_auto_fixture(n = 16L)

  for (type in c("fixed", "generalized_nn", "adaptive_nn")) {
    bw <- nplsqregbw(
      y ~ x1,
      data = dat,
      scale = rep.int(1, nrow(dat)),
      nomad = "auto",
      bwtype = type,
      degree.min = 0L,
      degree.max = 1L,
      nmulti = 1L,
      optim.control = list(maxit = 1L)
    )
    expect_identical(bw$reg.bws$type, type)
    expect_identical(bw$reg.bws$degree.search$mode, "exhaustive")
    expect_true(is.finite(bw$objective))
  }
})

test_that("nplsqreg AUTO resolves multidimensional smoothing to NOMAD/Powell", {
  options(np.messages = FALSE, np.tree = FALSE)
  dat <- make_nplsqreg_auto_fixture(p = 2L)
  bw <- nplsqregbw(
    y ~ x1 + x2,
    data = dat,
    scale = rep.int(1, nrow(dat)),
    nomad = "auto",
    degree.min = 0L,
    degree.max = 1L,
    nmulti = 1L,
    nomad.nmulti = 0L,
    optim.control = list(maxit = 1L),
    nomad.opts = list(MAX_BB_EVAL = 8L)
  )

  expect_identical(bw$reg.bws$degree.search$mode, "nomad+powell")
  expect_identical(bw$reg.bws$nomad.shortcut$source, "auto")
  expect_true(is.finite(bw$objective))
})

test_that("nplsqreg AUTO composes through vector tau and public fit methods", {
  options(np.messages = FALSE, np.tree = FALSE)
  dat <- make_nplsqreg_auto_fixture(n = 18L)
  fit <- nplsqreg(
    y ~ x1,
    data = dat,
    scale = rep.int(1, nrow(dat)),
    tau = c(0.4, 0.6),
    tau.search = "refined",
    nomad = "auto",
    degree.min = 0L,
    degree.max = 1L,
    nmulti = 1L,
    optim.control = list(maxit = 1L),
    residuals = TRUE
  )

  expect_s3_class(fit, "lsqregression")
  expect_true(all(vapply(
    fit$bws$tau.bws,
    function(x) identical(x$reg.bws$degree.search$mode, "exhaustive"),
    logical(1L)
  )))
  expect_identical(dim(fitted(fit)), c(nrow(dat), 2L))
  expect_identical(dim(residuals(fit)), c(nrow(dat), 2L))
  expect_silent(capture.output(summary(fit)))
  expect_silent(capture.output(summary(fit$bws)))
  expect_identical(dim(predict(fit, newdata = dat[1:3, "x1", drop = FALSE])),
                   c(3L, 2L))
})

test_that("nplsqreg AUTO rejects only incompatible controls", {
  dat <- make_nplsqreg_auto_fixture(n = 12L)
  common <- list(xdat = dat["x1"], ydat = dat$y,
                 scale = rep.int(1, nrow(dat)))

  expect_error(
    do.call(nplsqregbw, c(common, list(nomad = "invalid"))),
    "must be TRUE, FALSE, or \"auto\"",
    fixed = TRUE
  )
  expect_error(
    do.call(nplsqregbw, c(common, list(nomad = "auto", degree = 1L))),
    "does not support an explicit degree",
    fixed = TRUE
  )
  expect_error(
    do.call(nplsqregbw, c(common, list(nomad = "auto", nomad.pilot = "auto"))),
    "'nomad.pilot' must be TRUE or FALSE",
    fixed = TRUE
  )
  expect_error(
    do.call(nplsqregbw, c(common, list(
      nomad = "auto", degree.max = 1L, nomad.nmulti = 1L
    ))),
    "positive nomad.nmulti is only supported",
    fixed = TRUE
  )
  expect_error(
    do.call(nplsqregbw, c(common, list(
      nomad = "auto", bandwidth.compute = FALSE
    ))),
    "requires bandwidth.compute=TRUE",
    fixed = TRUE
  )
  expect_error(
    do.call(nplsqregbw, c(common, list(nomad = "auto", regtype = "lc"))),
    "requires regtype='lp'",
    fixed = TRUE
  )
  expect_error(
    do.call(nplsqregbw, c(common, list(nomad = "auto", basis = "tensor"))),
    "currently requires basis='glp'",
    fixed = TRUE
  )
  expect_error(
    do.call(nplsqregbw, c(common, list(
      nomad = TRUE, search.engine = "cell", degree.max = 1L
    ))),
    "requires search.engine='nomad' or 'nomad+powell'",
    fixed = TRUE
  )
})
