test_that("IV derivative initialization validates ties before preparation", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  z <- c(0, 0, .3, .6, 1, 1)
  w <- seq_along(z)/length(z)
  y <- z^2+w/10
  sentinel <- function(...) stop("preparation reached")
  # Probe the unchanged local body without broadcasting a master-only mock.
  owner <- getFromNamespace("npregivderiv.default","npRmpi")
  env <- new.env(parent=environment(owner))
  env$npudensbw <- sentinel
  env$.npRmpi_autodispatch_active <- function() FALSE
  environment(owner) <- env
  set.seed(4051)
  seed <- .Random.seed
  expect_error(owner(y=y,z=z,w=w,starting.values=seq_along(z)),
               "conflicting quadrature ordinates at duplicate coordinates")
  expect_identical(.Random.seed,seed)
  expect_error(owner(y=y,z=z,w=w,starting.values=2*z),
               "preparation reached")
  expect_error(owner(y=y,z=z,w=w), "preparation reached")
})
