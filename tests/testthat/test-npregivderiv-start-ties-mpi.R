test_that("conflicting IV starts fail coherently under public MPI dispatch", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages=FALSE)
  on.exit(options(old),add=TRUE)
  z <- c(0,0,.3,.6,1,1)
  w <- seq_along(z)/length(z)
  y <- z^2+w/10
  set.seed(4051)
  seed <- .Random.seed
  expect_error(npregivderiv(y=y,z=z,w=w,starting.values=seq_along(z)),
               "conflicting quadrature ordinates at duplicate coordinates")
  expect_identical(.Random.seed,seed)
  # An ordinary collective fit still works after the rejected request.
  fit <- npreg(txdat=data.frame(w=w),tydat=y,bws=.4,se=FALSE)
  expect_length(fitted(fit),length(y))
  expect_true(all(is.finite(fitted(fit))))
})
