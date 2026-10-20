test_that("distributed NN copula errors retain ratio and request contracts", {
  if (!spawn_mpi_slaves()) skip("Could not initialize MPI context")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  ns <- asNamespace("npRmpi")
  owner <- get(".npcopula_density_se_nn", ns)
  x <- data.frame(x = sin(1:24) + (1:24)/24,
                  y = cos(1:24) + (1:24)/48)
  ex <- data.frame(x = c(.15, .65, .15), y = c(.2, .6, .2))
  for (type in c("generalized_nn", "adaptive_nn")) {
    bw <- npudensbw(dat = x, bws = c(6, 6), bwtype = type,
                    bandwidth.compute = FALSE)
    answer <- owner(bw, x, ex)
    expect_true(all(is.finite(answer) & answer >= 0))
    expect_identical(answer[1L], answer[3L])
    expect_equal(owner(bw, x, ex[c(3, 2, 1), ]), answer[c(3, 2, 1)],
                  tolerance = 2e-10)
    single <- data.frame(x = ordered("a"))
    single.bw <- npudensbw(dat = single, bws = 0, bwtype = type,
                           bandwidth.compute = FALSE)
    expect_identical(owner(single.bw, single, single), 0)
    b1 <- npudensbw(dat = x[1L], bws = 6, bwtype = type,
                    bandwidth.compute = FALSE)
    expect_identical(owner(b1, x[1L], ex[1L]), rep(0, 3L))
    off <- npcopula(bws = bw, data = x)
    on <- npcopula(bws = bw, data = x, se = TRUE)
    expect_null(off$copulaerr)
    expect_identical(as.data.frame(off), as.data.frame(on))
    expect_identical(off$bws, on$bws)
    expect_identical(predict(on, se.fit = TRUE)$se.fit, se(on))
    set.seed(917)
    seed <- .Random.seed
    invisible(owner(bw, x, ex))
    expect_identical(.Random.seed, seed)
  }
})

test_that("NN ratio norm protects representable large and small errors", {
  norm <- get(".npcopula_density_se_nn_norm", asNamespace("npRmpi"))
  expect_equal(norm(c(-1e200, 1e200)), sqrt(2) * 1e200)
  expect_equal(norm(c(-1e-200, 1e-200)), sqrt(2) * 1e-200)
  expect_gt(norm(c(-1e-200, 1e-200)), 0)
  expect_identical(norm(c(4, 4)), 0)
  expect_identical(norm(c(NA_real_, 1)), NA_real_)
})
