test_that("explicit plot normals validate their request before allocation", {
  expect_error(.np_plot_normal_request(-1, c("Mersenne-Twister", "Box-Muller"), 1), "seed")
  expect_error(.np_plot_normal_request(42, c("user-supplied", "Box-Muller"), 1), "built-in")
  expect_error(.np_plot_normal_request(42, c("Mersenne-Twister", "Box-Muller"), .5), "count")
  stream <- .np_plot_normal_stream(42, c("Mersenne-Twister", "Box-Muller"), 0)
  on.exit(.np_plot_normal_close(stream), add=TRUE)
  expect_null(stream$path)
  expect_identical(.np_plot_normal_read(stream, 0), numeric())
  expect_error(.np_plot_normal_read(stream, 1), "read")
})

test_that("seeded plot normals preserve the caller's Box-Muller cache", {
  old.kind <- RNGkind()
  old.seed <- if (exists(".Random.seed", .GlobalEnv, inherits=FALSE)) .Random.seed else NULL
  on.exit({do.call(RNGkind, as.list(old.kind));
    if (is.null(old.seed)) rm(list=".Random.seed", envir=.GlobalEnv)
    else assign(".Random.seed", old.seed, .GlobalEnv)}, add=TRUE)
  RNGkind("Mersenne-Twister", "Box-Muller")
  set.seed(42); expected <- rnorm(9)
  set.seed(891); invisible(rnorm(1)); seed <- .Random.seed; following <- rnorm(7)
  set.seed(891); invisible(rnorm(1))
  stream <- .np_plot_normal_stream(42, RNGkind(), 9)
  on.exit(.np_plot_normal_close(stream), add=TRUE)
  expect_false(stream$info$pid == Sys.getpid())
  expect_false(any(c("np", "npRmpi", "Rmpi") %in% stream$info$namespaces))
  expect_identical(c(.np_plot_normal_read(stream, 1), .np_plot_normal_read(stream, 8)), expected)
  expect_identical(.Random.seed, seed)
  expect_identical(rnorm(7), following)
  path <- stream$path
  .np_plot_normal_close(stream)
  expect_false(dir.exists(path))
  expect_silent(.np_plot_normal_close(stream))
})
