test_that("bounded continuous copula margins retain physical bandwidths", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  d <- data.frame(x = seq(.05, .95, length.out = 24),
                  z = 1 + .4*sin(seq_len(24)))
  for (kernel in c("gaussian", "beta"))
    for (constructor in list(npudensbw, npudistbw)) {
      args <- list(dat = d, bws = c(.8, .9), bwscaling = TRUE,
        bandwidth.compute = FALSE, ckertype = kernel,
        ckerbound = "fixed", ckerlb = c(0, 0), ckerub = c(1, 2))
      b <- do.call(constructor, args)
      args$bws <- b$bandwidth$x
      args$bwscaling <- FALSE
      p <- do.call(constructor, args)
      a <- npcopula(b, data = d, se = TRUE)
      expected <- npcopula(p, data = d, se = TRUE)
      expect_identical(as.data.frame(a), as.data.frame(expected))
      expect_identical(se(a), se(expected))
    }
})
