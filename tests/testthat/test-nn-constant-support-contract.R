test_that("NN bandwidth owners reject constant continuous support before search", {
  pool.was.active <- .mpi_pool_active()
  skip_if_not(spawn_mpi_slaves(1L))
  on.exit(if (!pool.was.active) close_mpi_slaves(), add = TRUE)
  dat <- data.frame(x = rep(1, 8L), y = sin(seq_len(8L)))
  for (bt in c("generalized_nn", "adaptive_nn")) {
    expect_error(npudensbw(~x, data = dat, bwtype = bt, nmulti = 1L),
      "require at least two distinct continuous regressor values", fixed = TRUE)
    expect_error(npudistbw(~x, data = dat, bwtype = bt, nmulti = 1L),
      "require at least two distinct continuous regressor values", fixed = TRUE)
    expect_error(npregbw(y ~ x, data = dat, bwtype = bt, regtype = "lc", nmulti = 1L),
      "require at least two distinct continuous regressor values", fixed = TRUE)
    for (fun in list(npcdensbw, npcdistbw)) {
      expect_error(fun(y ~ x, data = dat, bwtype = bt, regtype = "lc", nmulti = 1L),
        "require at least two distinct continuous regressor values", fixed = TRUE)
      expect_error(fun(x ~ y, data = dat, bwtype = bt, regtype = "lc", nmulti = 1L),
        "require at least two distinct continuous variable values", fixed = TRUE)
    }
  }
})
