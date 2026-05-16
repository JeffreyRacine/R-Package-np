.conditional_fixed_formula_env <- new.env(parent = emptyenv())

.ensure_conditional_fixed_formula_pool <- function() {
  if (!isTRUE(.conditional_fixed_formula_env$started)) {
    npRmpi.init(nslaves = 1L, quiet = TRUE)
    .conditional_fixed_formula_env$started <- TRUE
    withr::defer({
      if (isTRUE(.conditional_fixed_formula_env$started)) {
        try(npRmpi.quit(force = TRUE), silent = TRUE)
        .conditional_fixed_formula_env$started <- FALSE
      }
    }, envir = testthat::teardown_env())
  }
}

test_that("conditional fixed formula routes honor numeric bws before autodispatch", {
  skip_on_cran()
  .ensure_conditional_fixed_formula_pool()
  old_opts <- options(np.messages = FALSE, np.tree = FALSE,
                      np.categorical.compress = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260623L)
  n <- 160L
  dat <- data.frame(
    y = ordered(sample(1:4, n, TRUE)),
    x = ordered(sample(1:4, n, TRUE))
  )

  bw.dens <- npcdensbw(y ~ x, data = dat, bws = c(.25, .25),
                       bandwidth.compute = FALSE)
  fit.dens.ref <- npcdens(bws = bw.dens)
  fit.dens.formula <- npcdens(y ~ x, data = dat, bws = c(.25, .25))

  bw.dist <- npcdistbw(y ~ x, data = dat, bws = c(.25, .25),
                       bandwidth.compute = FALSE)
  fit.dist.ref <- npcdist(bws = bw.dist)
  fit.dist.formula <- npcdist(y ~ x, data = dat, bws = c(.25, .25))

  expect_equal(fitted(fit.dens.formula), fitted(fit.dens.ref),
               tolerance = 1e-8)
  expect_equal(fitted(fit.dist.formula), fitted(fit.dist.ref),
               tolerance = 1e-8)
})

test_that("unconditional distribution formula route is preserved before autodispatch", {
  skip_on_cran()
  .ensure_conditional_fixed_formula_pool()
  old_opts <- options(np.messages = FALSE, np.tree = FALSE,
                      np.categorical.compress = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260624L)
  dat <- data.frame(x = ordered(sample(1:4, 80L, TRUE)))

  bw <- npudistbw(~ x, data = dat, nmulti = 1)
  fit.ref <- npudist(bws = bw)
  fit.formula <- npudist(~ x, data = dat, nmulti = 1)

  expect_equal(fitted(fit.formula), fitted(fit.ref), tolerance = 1e-8)
})
