test_that(".np_seed_exit restores prior seed when one exists", {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    seed_orig <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", seed_orig, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(123)
  seed_before <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  state <- .np_seed_enter(42)
  .np_seed_exit(state)
  seed_after <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  expect_identical(seed_after, seed_before)
})

test_that(".np_seed_exit remove_if_absent drops synthetic seed", {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    seed_orig <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    rm(".Random.seed", envir = .GlobalEnv)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", seed_orig, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  state <- .np_seed_enter(42)
  expect_true(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  .np_seed_exit(state, remove_if_absent = TRUE)
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})
