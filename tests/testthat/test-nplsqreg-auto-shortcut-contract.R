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
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  dat <- make_nplsqreg_auto_fixture(n = 1024L)
  npseed(20260813L)
  bw <- nplsqregbw(
    y ~ x1,
    data = dat,
    scale = rep.int(1, nrow(dat)),
    nomad = "auto",
    degree.min = 0L,
    degree.max = 1L,
    nmulti = 1L,
    optim.control = list(maxit = 1L)
  )

  expect_s3_class(bw, "lsqregressionbandwidth")
  expect_identical(bw$reg.bws$degree.search$mode, "exhaustive")
  expect_true(isTRUE(bw$reg.bws$degree.search$certified))
  expect_identical(bw$reg.bws$nomad.shortcut$source, "auto")
  expect_identical(bw$reg.bws$nomad.shortcut$nomad, "auto")
  expect_true(is.finite(bw$delta))
  expect_true(is.finite(bw$objective))
})

test_that("nplsqreg AUTO isolates every MPI cell in a public transaction", {
  source <- paste(deparse(getFromNamespace(".nplsqreg_cell_search", "npRmpi")),
                  collapse = "\n")
  fixed <- paste(deparse(getFromNamespace(".npRmpi_nplsqreg_fixed_cell", "npRmpi")),
                 collapse = "\n")
  gate <- paste(deparse(getFromNamespace(
    ".npRmpi_nplsqreg_should_coordinate_cell", "npRmpi"
  )), collapse = "\n")
  gate.compact <- gsub("[[:space:]]+", "", gate)

  expect_match(source, ".npRmpi_nplsqreg_fixed_cell", fixed = TRUE)
  expect_match(source, ".nplsqreg_call_fixed_degree_core", fixed = TRUE)
  expect_match(source, "bandwidth.compute = TRUE", fixed = TRUE)
  expect_match(source, ".np_degree_search", fixed = TRUE)
  expect_match(source, "direction = \"min\"", fixed = TRUE)
  expect_match(fixed, "do.call(nplsqregbw, fixed.args)", fixed = TRUE)
  expect_match(fixed, "nomad = FALSE", fixed = TRUE)
  expect_match(gate, "identical(nomad, \"auto\")", fixed = TRUE)
  expect_match(gate.compact, "identical(sum(xtype$icon),1L)", fixed = TRUE)
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
