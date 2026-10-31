test_that("regression percentile tasks replay the coordinator generator", {
  kind <- RNGkind(); had <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  old <- if (had) get(".Random.seed", .GlobalEnv) else NULL
  on.exit({do.call(RNGkind, as.list(kind)); if (had) assign(".Random.seed", old, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)}, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  state <- .np_reg_percentile_rng_state(91)
  contract <- list(runtime = paste(R.version$major, R.version$minor, sep = "."),
    header = state[1L], length = length(state), kind = RNGkind())
  expected <- .np_reg_percentile_with_seed(91, runif(37))
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection"); set.seed(17)
  worker <- .Random.seed
  expect_identical(.np_reg_percentile_with_state(state, contract, runif(37)), expected)
  expect_identical(.Random.seed, worker)
  expect_error(.np_reg_percentile_with_state(NULL, contract, stop("entered")), "invalid regression percentile")
  expect_error(.np_reg_percentile_with_state(state, contract, stop("primary")), "primary")
  expect_identical(.Random.seed, worker)
  state[20L] <- NA_integer_
  expect_silent(.np_reg_percentile_rng_validate(state, contract))
})
