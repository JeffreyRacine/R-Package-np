test_that("adaptive density CV rejects only out-of-domain fold candidates", {
  skip_on_cran()
  if (!spawn_mpi_slaves(1L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE, np.tree = FALSE, np.extendednn = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = (1:9)^2 / 81)
  for (method in c("cv.ml", "cv.ls")) {
    for (k in 6:8) {
      bw <- npudensbw(dat = x, bws = k, bwtype = "adaptive_nn",
        bwmethod = method, bandwidth.compute = FALSE)
      raw <- npRmpi:::npudensbw.bandwidth(dat = x, bws = bw,
        eval.only = TRUE, invalid.penalty = "dbmax", nmulti = 1L)
      if (k == 8L) expect_identical(as.double(raw$fval), -.Machine$double.xmax)
      else expect_true(is.finite(raw$fval) && abs(raw$fval) < .Machine$double.xmax)
      # Full-sample fitting still has eight donor neighbours available.
      fit <- npudens(bws = bw, tdat = x)
      expect_true(all(is.finite(fitted(fit))))
    }
  }
})

test_that("adaptive density Powell survives a rejected fold boundary trial", {
  skip_on_cran()
  if (!spawn_mpi_slaves(1L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE, np.tree = FALSE, np.extendednn = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = c(rep(0, 20L), 1:4))
  set.seed(9203L)
  bw <- npudensbw(dat = x, bwtype = "adaptive_nn", bwsolver = "powell",
    nmulti = 1L, itmax = 20L, powell.remin = FALSE)
  raw <- npRmpi:::npudensbw.bandwidth(dat = x, bws = bw,
    eval.only = TRUE, invalid.penalty = "dbmax", nmulti = 1L)
  expect_true(is.finite(bw$fval) && abs(bw$fval) < .Machine$double.xmax)
  expect_equal(as.double(bw$fval), as.double(raw$fval), tolerance = 2e-12)
})
