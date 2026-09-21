test_that("smooth coefficient positional roles follow the documented signature", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force=TRUE), add=TRUE)
  old <- options(np.messages=FALSE)
  on.exit(options(old), add=TRUE)
  x <- seq(-1, 1, length.out=24)
  y <- sin(2*x)+cos(seq_along(x))/10
  z <- x^3
  named <- npscoef(bws=.45, txdat=x, tydat=y, se=FALSE)
  positional <- npscoef(.45, x, y, se=FALSE)
  expect_equal(positional$bws$bw, .45, tolerance=0)
  expect_identical(fitted(positional), fitted(named))
  set.seed(730)
  direct <- npscoef(txdat=x, tydat=y, tzdat=z, nmulti=1,
                    optim.maxit=10, se=FALSE)
  rng <- .Random.seed
  set.seed(730)
  data.first <- npscoef(x, y, zdat=z, nmulti=1,
                        optim.maxit=10, se=FALSE)
  expect_identical(.Random.seed, rng)
  expect_identical(data.first$bws$bw, direct$bws$bw)
  expect_identical(fitted(data.first), fitted(direct))
  a <- npscoefbw(x, y, z, bws=.45, bandwidth.compute=FALSE)
  b <- npscoefbw(xdat=x, ydat=y, zdat=z, bws=.45, bandwidth.compute=FALSE)
  expect_equal(a$bw, .45, tolerance=0)
  expect_identical(a$bw, b$bw)
})
