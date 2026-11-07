# Small pure contract: broad custom-generator/MPI matrices live outside default lane.
test_that("npsig replay checks only assigned uniform generator headers", {
  ns <- asNamespace("npRmpi")
  owner <- get("npsigtest.rbandwidth", ns)
  check <- get(".npRmpi_npsig_replay_check", environment(owner))
  before <- if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
    get(".Random.seed", .GlobalEnv) else NULL
  plans <- list(c(10403L, 1L), c(10203L, 1L), c(10303L, 1L),
                c(10405L, 1L), c(10305L, 1L))
  expect_invisible(check(1:3, plans))
  expect_invisible(check(integer(0), plans))
  expect_error(check(c(1L, 4L), plans), "user-supplied generator state")
  expect_error(check(5L, plans), "user-supplied generator state")
  after <- if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
    get(".Random.seed", .GlobalEnv) else NULL
  expect_identical(after, before)
})
