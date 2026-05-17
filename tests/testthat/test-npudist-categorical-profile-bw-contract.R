test_that("npudist ordered-categorical profile bandwidth CV matches dense CV", {
  old.tree <- getOption("np.tree")
  old.compress <- getOption("np.categorical.compress")
  old.messages <- getOption("np.messages")
  on.exit({
    options(np.tree = old.tree)
    options(np.categorical.compress = old.compress)
    options(np.messages = old.messages)
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

    expect_equal(profile$bw, dense$bw, tolerance = 1e-12)
    expect_equal(profile$fval, dense$fval, tolerance = 1e-12)
    expect_equal(profile$num.feval, dense$num.feval)
  }

  run_case("wangvanryzin", 20260538L)
  run_case("liracine", 20260539L)
  run_case("racineliyan", 20260540L)
})
