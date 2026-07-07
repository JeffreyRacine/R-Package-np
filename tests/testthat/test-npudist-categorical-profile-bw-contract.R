test_that("npudist ordered-categorical profile bandwidth CV matches dense CV", {
  if (exists("spawn_mpi_slaves", mode = "function"))
    spawn_mpi_slaves(1L)

  old.tree <- getOption("np.tree")
  old.compress <- getOption("np.categorical.compress")
  old.messages <- getOption("np.messages")
  on.exit({
    options(np.tree = old.tree)
    options(np.categorical.compress = old.compress)
    options(np.messages = old.messages)
    if (exists("close_mpi_slaves", mode = "function"))
      close_mpi_slaves()
  }, add = TRUE)
  options(np.messages = FALSE)

  run_case <- function(okertype, seed) {
    set.seed(seed)
    n <- 240L
    dat <- data.frame(
      o1 = ordered(sample(1:6, n, TRUE)),
      o2 = ordered(sample(1:5, n, TRUE)),
      o3 = ordered(sample(1:4, n, TRUE))
    )

    options(np.tree = FALSE, np.categorical.compress = FALSE)
    dense <- npudistbw(
      ~ o1 + o2 + o3,
      data = dat,
      nmulti = 1,
      okertype = okertype
    )

    options(np.tree = FALSE, np.categorical.compress = TRUE)
    profile <- npudistbw(
      ~ o1 + o2 + o3,
      data = dat,
      nmulti = 1,
      okertype = okertype
    )

    expect_equal(profile$fval, dense$fval, tolerance = 1e-8)
    expect_gt(profile$num.feval.fast, 0)
  }

  run_case("wangvanryzin", 20260538L)
  run_case("liracine", 20260539L)
  run_case("racineliyan", 20260540L)
})

test_that("npudist ordered-categorical profile bandwidth CV supports training-grid integral", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(), add = TRUE)

  old.tree <- getOption("np.tree")
  old.compress <- getOption("np.categorical.compress")
  old.messages <- getOption("np.messages")
  on.exit({
    options(np.tree = old.tree)
    options(np.categorical.compress = old.compress)
    options(np.messages = old.messages)
  }, add = TRUE)
  options(np.messages = FALSE)

  set.seed(20260541L)
  n <- 360L
  dat <- data.frame(
    o1 = ordered(sample(1:6, n, TRUE)),
    o2 = ordered(sample(1:5, n, TRUE))
  )

  options(np.tree = FALSE, np.categorical.compress = FALSE)
  dense <- npudistbw(
    ~ o1 + o2,
    data = dat,
    nmulti = 1,
    do.full.integral = TRUE
  )

  options(np.tree = FALSE, np.categorical.compress = TRUE)
  profile <- npudistbw(
    ~ o1 + o2,
    data = dat,
    nmulti = 1,
    do.full.integral = TRUE
  )

  expect_equal(profile$fval, dense$fval, tolerance = 1e-8)
  expect_gt(profile$num.feval.fast, 0)
})

test_that("npudist generated ordered grids preserve interval-labelled levels", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(), add = TRUE)

  old.tree <- getOption("np.tree")
  old.compress <- getOption("np.categorical.compress")
  old.messages <- getOption("np.messages")
  on.exit({
    options(np.tree = old.tree)
    options(np.categorical.compress = old.compress)
    options(np.messages = old.messages)
  }, add = TRUE)
  options(np.messages = FALSE, np.tree = FALSE,
          np.categorical.compress = TRUE)

  set.seed(20260542L)
  n <- 180L
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  dat <- data.frame(
    o1 = ordered(cut(x1, quantile(x1, seq(0, 1, length.out = 5)),
                     include.lowest = TRUE)),
    o2 = ordered(cut(x2, quantile(x2, seq(0, 1, length.out = 5)),
                     include.lowest = TRUE))
  )

  expect_warning(
    bw <- npudistbw(~ o1 + o2, data = dat, nmulti = 1,
                    okertype = "liracine"),
    NA
  )
  expect_true(is.finite(bw$fval))
})
