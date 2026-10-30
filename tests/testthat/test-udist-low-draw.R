test_that("defined low-draw CDF bands are returned with one resolution warning", {
  old <- options(np.messages = FALSE, np.plot.progress = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = seq(-2.5, 2.5, length.out = 31L))
  bw <- npudistbw(dat = x, bws = .45, bandwidth.compute = FALSE)
  for (law in c("inid", "fixed", "geom")) {
    conditions <- list()
    z <- withCallingHandlers(plot(bw, xdat = x, errors = "bootstrap",
      bootstrap = law, center = "bias-corrected", band = "all", B = 1L,
      boot.control = np_boot_control(blocklen = if (law == "inid") NULL else 3L),
      neval = 5L, output = "data", random.seed = 99102905),
      np_band_resolution_warning = function(w) {
        conditions[[length(conditions) + 1L]] <<- w
        invokeRestart("muffleWarning")
      })$d1$inference
    expect_length(conditions, 1L)
    expect_named(z$bands, c("pointwise", "simultaneous", "bonferroni"))
    expect_true(all(is.finite(unlist(z$bands))))
    expect_true(all(z$resolution$computed))
    expect_true(all(z$resolution$resolution_limited))
  }
})

test_that("one numeric GNN draw is an undefined empirical SD, not a resolution refusal", {
  x <- data.frame(x = seq(-2.5, 2.5, length.out = 31L))
  bw <- npudistbw(dat = x, bws = 10L, bwtype = "generalized_nn",
                  bandwidth.compute = FALSE)
  state <- .np_udist_bc_prepare(bw, x)
  set.seed(99102902); before <- .Random.seed
  expect_error(.np_udist_joint_inference(state, x[1:3, , drop = FALSE],
    law = "inid", B = 1L), "empirical standard deviation requires at least two")
  expect_identical(.Random.seed, before)
  expect_error(.np_udist_joint_inference(state, x[1:3, , drop = FALSE],
    law = "inid", B = 0L), "CDF draws")
  expect_identical(.Random.seed, before)
})
