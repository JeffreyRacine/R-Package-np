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
    categorical <- data.frame(x = ordered(rep(letters[1:3], 4L)))
    categorical.bw <- npudensbw(dat = categorical, bws = .2, bwtype = type,
      ckertype = "uniform", bandwidth.compute = FALSE)
    expect_silent(categorical.se <- owner(categorical.bw, categorical, categorical))
    expect_identical(categorical.se, rep(0, nrow(categorical)))

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

test_that("retained NN metadata preserves genuine kernel conditions", {
  ns <- asNamespace(getNamespaceName(environment(npcopula)))
  retained <- get(".npcopula_density_se_nn_retained", ns)
  expect_warning(retained(list(ckertype = "uniform"), warning("distinct condition")),
    "distinct condition")
  expect_error(retained(list(ckertype = "uniform"), stop("distinct failure")),
    "distinct failure")
  advisory <- unname(get(".np_io_prefix_text", ns)(
    "ignoring kernel order specified with uniform kernel type"))
  expect_warning(retained(list(ckertype = "gaussian"), warning(advisory)),
    "ignoring kernel order")
})


test_that("uniform ANN query factors reuse invocation-local training preparation", {
  if (exists("spawn_mpi_slaves", mode = "function")) {
    if (!spawn_mpi_slaves()) skip("Could not initialize MPI context")
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }
  ns <- asNamespace(getNamespaceName(environment(npcopula)))
  x <- data.frame(x = sin(1:32) + (1:32)/32,
                  y = cos(1:32) + (1:32)/64)
  bw <- npudensbw(dat = x, bws = c(6, 6), bwtype = "adaptive_nn",
                  ckertype = "uniform", bandwidth.compute = FALSE)
  retained <- get(".npcopula_density_se_nn_retained", ns)
  prepare <- get(".np_plot_density_ann_prepare", ns)
  factor <- get(".npcopula_density_se_nn_factor", ns)
  context <- retained(bw, prepare(bw, x, TRUE, uniform.faces = TRUE))
  before <- serialize(context, NULL)
  evaluation <- rbind(x[c(1, 5, 11), ], data.frame(x = 1e6, y = 1e6))
  for (i in seq_len(nrow(evaluation)))
    expect_identical(factor(bw, x, evaluation[i, , drop = FALSE], context),
                     factor(bw, x, evaluation[i, , drop = FALSE]))
  expect_identical(serialize(context, NULL), before)
  # Fresh preparation must continue to reject invalid training geometry.
  tied <- x; tied$x <- 1
  expect_error(prepare(bw, tied, TRUE), "nonregular ANN spacing geometry")
})
