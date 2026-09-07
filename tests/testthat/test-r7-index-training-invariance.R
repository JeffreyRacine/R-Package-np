test_that("GNN index training statistics are independent of evaluation requests", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(11)
  n <- 60L
  x <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  signal <- with(x, x1 + .5 * x2 - .3 * x3)
  noise <- rnorm(n, sd = .8)
  for (method in c("ichimura", "kleinspady")) {
    y <- if (method == "ichimura") signal + noise else
      as.double(signal + noise > 0)
    for (regtype in c("lc", "ll", "lp")) {
      bw.args <- list(xdat = x, ydat = y, method = method,
        regtype = regtype,
        bwtype = "generalized_nn", bws = c(1, .5, -.3, 15),
        bandwidth.compute = FALSE)
      if (regtype == "lp") bw.args$degree <- 2L
      bw <- do.call(npindexbw, bw.args)
      ref <- npindex(bws = bw, txdat = x, tydat = y, residuals = TRUE)
      for (mode in list(list(), list(gradients = TRUE),
                        list(se.type = "bootstrap", B = 3L),
                        list(se = FALSE, gradients = TRUE))) {
        set.seed(123)
        a <- do.call(npindex, c(list(bws = bw, txdat = x, tydat = y,
                                    residuals = TRUE), mode))
        set.seed(123)
        b <- do.call(npindex, c(list(bws = bw, txdat = x, tydat = y,
                                    exdat = x[1:5, ], residuals = TRUE), mode))
        expect_equal(a$R2, b$R2, tolerance = 1e-12)
        expect_equal(a$resid, b$resid, tolerance = 1e-12)
        expect_equal(a$fit.mcfadden, b$fit.mcfadden, tolerance = 1e-12)
        expect_equal(fitted(b), fitted(a)[1:5], tolerance = 1e-12)
        if (!identical(mode$se, FALSE)) {
          expect_equal(vcov(a), vcov(b), tolerance = 1e-12)
          expect_equal(vcov(a), vcov(ref), tolerance = 1e-12)
        }
      }
    }
  }
})
