test_that("regression profile owners propagate original kernel errors", {
  if (exists("spawn_mpi_slaves", mode = "function")) {
    if (!spawn_mpi_slaves()) skip("Could not initialize MPI context")
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }
  pkg <- getNamespaceName(environment(npregbw))
  old <- options(np.messages = FALSE, np.categorical.compress = TRUE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = factor(rep(c("a", "b"), 20L)))
  y <- sin(seq_len(40L))
  bw <- npregbw(xdat = x, ydat = y, bws = .2,
                regtype = "lc", bandwidth.compute = FALSE)
  mean.owner <- getFromNamespace(".np_regression_cat_profile_mean", pkg)
  boot.owner <- getFromNamespace(".np_regression_cat_profile_boot_setup", pkg)
  kernel <- getFromNamespace(".np_regression_cat_profile_kernel_matrix", pkg)
  state <- new.env(parent = emptyenv())
  state$calls <- 0L
  state$fail.at <- 1L
  witness <- structure(list(message = "profile kernel witness", call = NULL),
                       class = c("profile_witness", "error", "condition"))
  local_mocked_bindings(.np_regression_cat_profile_kernel_matrix = function(...) {
    state$calls <- state$calls + 1L
    if (state$calls == state$fail.at) stop(witness)
    kernel(...)
  }, .package = pkg)
  expect_error(mean.owner(bw, x, y, x[1:2, , drop = FALSE]),
               "profile kernel witness", class = "profile_witness")
  for (at in 1:2) {
    state$calls <- 0L
    state$fail.at <- at
    expect_error(boot.owner(x, x[1:2, , drop = FALSE], y, bw),
                 "profile kernel witness", class = "profile_witness")
    expect_identical(state$calls, at)
  }
  state$calls <- 0L
  state$fail.at <- 1L
  state$alternate <- 0L
  local_mocked_bindings(.np_inid_boot_from_regression = function(...) {
    state$alternate <- state$alternate + 1L
    stop("alternate owner executed")
  }, .package = pkg)
  expect_error(plot(bw, output = "data", errors = "bootstrap",
                    bootstrap = "inid", B = 39L, neval = 7L),
               "profile kernel witness", class = "profile_witness")
  expect_identical(state$alternate, 0L)
})

test_that("constant-basis bootstrap preserves the original hat error", {
  if (exists("spawn_mpi_slaves", mode = "function")) {
    if (!spawn_mpi_slaves()) skip("Could not initialize MPI context")
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }
  pkg <- getNamespaceName(environment(npregbw))
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = seq(0, 1, length.out = 40L))
  y <- sin(seq_len(40L))
  bw <- npregbw(xdat = x, ydat = y, bws = .2,
                regtype = "lc", bandwidth.compute = FALSE)
  owner <- getFromNamespace(".np_inid_boot_from_regression", pkg)
  state <- new.env(parent = emptyenv())
  state$alternate <- 0L
  witness <- structure(list(message = "hat construction witness", call = NULL),
                       class = c("hat_witness", "error", "condition"))
  local_mocked_bindings(npreghat.rbandwidth = function(...) stop(witness),
    .np_inid_boot_from_regression_localpoly_fixed = function(...) {
      state$alternate <- state$alternate + 1L
      stop("alternate owner executed")
    }, .package = pkg)
  expect_error(owner(xdat = x, exdat = x[1:7, , drop = FALSE], bws = bw,
                     ydat = y, B = 2L, counts = matrix(1, 40L, 2L)),
               "hat construction witness", class = "hat_witness")
  expect_identical(state$alternate, 0L)
})
