test_that("npudist all-ordered profile route preserves fitted distribution", {
  options(np.messages = FALSE)

  for (okertype in c("liracine", "racineliyan")) {
    set.seed(20260537)
    n <- 512
    dat <- data.frame(
      o1 = ordered(sample(1:6, n, TRUE)),
      o2 = ordered(sample(1:5, n, TRUE)),
      o3 = ordered(sample(1:4, n, TRUE))
    )

    old <- options(np.tree = FALSE, np.categorical.compress = FALSE)
    dense_bw <- npudistbw(
      ~ o1 + o2 + o3,
      data = dat,
      bws = c(0.15, 0.25, 0.35),
      bandwidth.compute = FALSE,
      okertype = okertype
    )
    dense <- npudist(bws = dense_bw)
    options(old)

    old <- options(np.tree = FALSE, np.categorical.compress = TRUE)
    profile_bw <- npudistbw(
      ~ o1 + o2 + o3,
      data = dat,
      bws = c(0.15, 0.25, 0.35),
      bandwidth.compute = FALSE,
      okertype = okertype
    )
    profile <- npudist(bws = profile_bw)
    options(old)

    expect_equal(fitted(profile), fitted(dense), tolerance = 1e-12)
    expect_equal(se(profile), se(dense), tolerance = 1e-12)
  }
})
