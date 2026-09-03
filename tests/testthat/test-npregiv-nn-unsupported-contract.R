test_that("proof-of-concept IV estimators reject NN before computation", {
  message <- paste0("npregiv() and npregivderiv() do not currently support ",
                    "nearest-neighbour bandwidths; use bwtype = \"fixed\" (the default).")
  dat <- data.frame(y = seq_len(8L), z = seq_len(8L) / 9,
                    w = sin(seq_len(8L)))
  for (fun in list(npregiv, npregivderiv)) {
    for (bt in c("generalized_nn", "adaptive_nn", "g", "a")) {
      set.seed(923L)
      before <- .Random.seed
      expect_error(fun(y = dat$y, z = dat$z, w = dat$w,
                       bwtype = bt, ckertype = stop("kernel argument forced")),
                   message, fixed = TRUE)
      expect_identical(.Random.seed, before)
      expect_error(fun(y ~ z | w, data = dat, bwtype = bt), message, fixed = TRUE)
      expect_identical(.Random.seed, before)
    }
  }
})

test_that("the IV bandwidth guard preserves lazy unrelated arguments", {
  guard <- get(".np_iv_require_fixed_bandwidth", envir = environment(npregiv))
  expect_null(guard())
  expect_null(guard(bwtype = "fixed"))
  expect_null(guard(bwtype = "f"))
  expect_null(guard(ckertype = stop("forced")))
  expect_null(guard(stop("unnamed argument forced")))
})
