transform_bounds_fixture <- function() {
  n <- 36L
  x <- seq(0.03, 0.97, length.out = n)
  u <- factor(rep(c("a", "b", "c"), length.out = n))
  o <- ordered(rep(1:4, length.out = n))
  y <- sin(2 * pi * x) + c(a = 0, b = 0.25, c = -0.15)[u] +
    0.08 * (as.integer(o) - 2.5)
  data.frame(y, x, u, o)
}

transform_search_snapshot <- function(fun, args, seed = 314159L) {
  npseed(seed)
  bw <- do.call(fun, args)
  list(
    fval = bw[["fval", exact = TRUE]],
    fval.history = bw[["fval.history", exact = TRUE]],
    eval.history = bw[["eval.history", exact = TRUE]],
    invalid.history = bw[["invalid.history", exact = TRUE]],
    bw = bw[["bw", exact = TRUE]],
    xbw = bw[["xbw", exact = TRUE]],
    ybw = bw[["ybw", exact = TRUE]],
    sfactor = bw[["sfactor", exact = TRUE]]
  )
}

test_that("regression transformed search starts on the external bandwidth scale", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  dat <- transform_bounds_fixture()
  set.seed(314159)
  bw <- npregbw(
    y ~ x + u + o, data = dat,
    scale.factor.init = 0.55, nmulti = 1L,
    transform.bounds = TRUE, powell.remin = FALSE, itmax = 1L
  )

  expect_true(all(is.finite(bw$bw)))
  expect_true(all(bw$bw > 0))
  expect_true(is.finite(bw$fval))

  replay <- npRmpi:::.npregbw_eval_only(dat[c("x", "u", "o")], dat$y, bw)
  expect_equal(replay$objective, bw$fval, tolerance = 64 * .Machine$double.eps)
})

test_that("conditional-distribution later starts use transformed coordinates", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  dat <- transform_bounds_fixture()
  set.seed(271828)
  bw <- npcdistbw(
    y ~ x + u, data = dat,
    scale.factor.init = 0.55, nmulti = 2L,
    transform.bounds = TRUE, powell.remin = FALSE, itmax = 1L
  )

  expect_true(all(is.finite(unlist(bw$sfactor))))
  expect_true(all(unlist(bw$sfactor) > 0))
  expect_true(is.finite(bw$fval))
  expect_length(bw$fval.history, 2L)
  expect_true(all(is.finite(bw$fval.history)))
  expect_equal(bw$fval, min(bw$fval.history), tolerance = 1e-12)

  replay <- npRmpi:::.npcdistbw_eval_only(dat[c("x", "u")], dat$y, bws = bw)
  expect_equal(replay$objective, bw$fval, tolerance = 5e-7)
})

test_that("transformed first starts do not inherit native search state", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  dat <- transform_bounds_fixture()
  common <- list(
    scale.factor.init = 0.55, nmulti = 1L,
    transform.bounds = TRUE, powell.remin = FALSE, itmax = 1L
  )
  owners <- list(
    npudens = list(
      fun = npudensbw,
      args = c(list(dat = dat[c("x", "u", "o")]), common)
    ),
    npudist = list(
      fun = npudistbw,
      args = c(list(dat = dat["x"]), common)
    ),
    npcdens = list(
      fun = npcdensbw,
      args = c(list(formula = y ~ x + u + o, data = dat), common)
    ),
    npcdist = list(
      fun = npcdistbw,
      args = c(list(formula = y ~ x + u + o, data = dat), common)
    )
  )

  .Call("C_np_release_static_buffers", PACKAGE = "npRmpi")
  for (owner in owners) {
    first <- transform_search_snapshot(owner$fun, owner$args)
    second <- transform_search_snapshot(owner$fun, owner$args)
    expect_identical(second, first)
  }

  .Call("C_np_release_static_buffers", PACKAGE = "npRmpi")
  target <- transform_search_snapshot(owners$npcdist$fun, owners$npcdist$args)
  invisible(transform_search_snapshot(owners$npudens$fun, owners$npudens$args))
  after_other_owner <- transform_search_snapshot(
    owners$npcdist$fun, owners$npcdist$args
  )
  expect_identical(after_other_owner, target)
})

test_that("conditional transformed starts respect split unordered kernels", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  dat <- transform_bounds_fixture()
  dat$z <- factor(rep(c("low", "high"), each = 3L, length.out = nrow(dat)))
  args <- list(
    formula = z ~ u + x, data = dat,
    uykertype = "aitchisonaitken", uxkertype = "liracine",
    scale.factor.init = 0.55, nmulti = 1L,
    transform.bounds = TRUE, powell.remin = FALSE, itmax = 1L
  )

  .Call("C_np_release_static_buffers", PACKAGE = "npRmpi")
  first <- transform_search_snapshot(npcdensbw, args)
  second <- transform_search_snapshot(npcdensbw, args)
  expect_identical(second, first)
})
