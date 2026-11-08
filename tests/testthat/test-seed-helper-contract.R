test_that("temporary scopes restore supported caller state on return and error", {
  set.seed(123)
  before <- .Random.seed
  expected <- c(runif(4), rnorm(5))
  assign(".Random.seed", before, .GlobalEnv)
  expect_identical(.np_with_seed(42L, invisible(2L)), 2L)
  expect_error(.np_with_seed(42L, {rnorm(3); stop("body error")}), "body error")
  expect_identical(.Random.seed, before)
  expect_identical(c(runif(4), rnorm(5)), expected)
})

test_that("an absent supported caller seed is restored on success and error", {
  kind <- RNGkind(); set.seed(123); before <- .Random.seed
  on.exit({do.call(RNGkind, as.list(kind)); assign(".Random.seed", before, .GlobalEnv)}, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion"); set.seed(42L)
  expected <- c(runif(3), rnorm(2))
  rm(".Random.seed", envir = .GlobalEnv)
  expect_identical(.np_with_seed(42L, c(runif(3), rnorm(2))), expected)
  expect_false(exists(".Random.seed", .GlobalEnv, inherits = FALSE))
  expect_error(.np_with_seed(42L, stop("body error")), "body error")
  expect_false(exists(".Random.seed", .GlobalEnv, inherits = FALSE))
})

test_that("cached normals are refused without changing their continuation", {
  kind <- RNGkind(); set.seed(123); before <- .Random.seed
  on.exit({do.call(RNGkind, as.list(kind)); assign(".Random.seed", before, .GlobalEnv)}, add=TRUE)
  RNGkind(normal.kind = "Box-Muller"); set.seed(21L)
  invisible(rnorm(1)); expected <- rnorm(1)
  set.seed(21L); invisible(rnorm(1)); visible <- .Random.seed
  expect_error(.np_with_seed(42L, stop("body entered")), "Box-Muller")
  expect_identical(.Random.seed, visible)
  expect_identical(rnorm(1), expected)
  set.seed(21L); invisible(rnorm(1)); visible <- .Random.seed
  rm(".Random.seed", envir = .GlobalEnv)
  expect_error(.np_with_seed(42L, stop("body entered")), "Box-Muller")
  expect_error(.np_rng_scope_capture(uniform.only = TRUE), "Box-Muller")
  expect_false(exists(".Random.seed", .GlobalEnv, inherits = FALSE))
  assign(".Random.seed", visible, .GlobalEnv)
  expect_identical(rnorm(1), expected)
})

test_that("uniform-only transport does not certify a generator reset", {
  set.seed(123); before <- .Random.seed
  scope <- .np_rng_scope_capture(uniform.only = TRUE)
  expect_error(.np_rng_scope_seed(scope, 42L), "uniform-only")
  expect_error(.np_rng_scope_state(scope, before), "uniform-only")
  expect_identical(.Random.seed, before)
})
