library(np)

make_jksum_mixed_data <- function(n = 220L, seed = 42L) {
  set.seed(seed)
  x1 <- runif(n)
  x2 <- runif(n)
  z1 <- rbinom(n, 1, 0.5)
  z2 <- rbinom(n, 1, 0.5)
  y <- cos(2 * pi * x1) + 0.5 * sin(2 * pi * x2) + z1 + rnorm(n, sd = 0.20)
  data.frame(
    y = y,
    x1 = x1,
    x2 = x2,
    z1 = factor(z1),
    z2 = ordered(z2)
  )
}

run_reg_cv_once <- function(dat, regtype, bwmethod) {
  t_bw <- system.time(
    bw <- npregbw(
      y ~ x1 + x2 + z1 + z2,
      regtype = regtype,
      bwmethod = bwmethod,
      nmulti = 1,
      data = dat
    )
  )

  list(
    fval = as.numeric(bw$fval),
    nfe = as.integer(bw$num.feval),
    elapsed = as.numeric(t_bw[["elapsed"]])
  )
}

test_that("jksum regression CV parity is deterministic for mixed data", {
  old_opts <- options(
    np.messages = FALSE,
    np.largeh.rel.tol = 0.05,
    np.disc.upper.rel.tol = 0.05
  )
  on.exit(options(old_opts), add = TRUE)

  dat <- make_jksum_mixed_data(n = 220L, seed = 100L)
  combos <- expand.grid(
    regtype = c("ll", "lc"),
    bwmethod = c("cv.ls", "cv.aic"),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(combos))) {
    regtype <- combos$regtype[[i]]
    bwmethod <- combos$bwmethod[[i]]

    set.seed(123)
    r1 <- run_reg_cv_once(dat, regtype, bwmethod)
    set.seed(123)
    r2 <- run_reg_cv_once(dat, regtype, bwmethod)

    expect_true(is.finite(r1$fval))
    expect_true(is.finite(r2$fval))
    expect_true(r1$nfe > 0L)
    expect_true(r2$nfe > 0L)
    expect_equal(r1$nfe, r2$nfe)
    expect_equal(r1$fval, r2$fval, tolerance = 1e-12)
  }
})

test_that("jksum regression CV smoke performance remains bounded", {
  skip_on_cran()

  old_opts <- options(
    np.messages = FALSE,
    np.largeh.rel.tol = 0.05,
    np.disc.upper.rel.tol = 0.05
  )
  on.exit(options(old_opts), add = TRUE)

  dat <- make_jksum_mixed_data(n = 220L, seed = 101L)

  set.seed(321)
  r_ll <- run_reg_cv_once(dat, "ll", "cv.ls")
  set.seed(321)
  r_lc <- run_reg_cv_once(dat, "lc", "cv.ls")

  expect_true(r_ll$nfe > 0L)
  expect_true(r_lc$nfe > 0L)
  expect_true(is.finite(r_ll$fval))
  expect_true(is.finite(r_lc$fval))

  # Guardrail: this should remain a small smoke test in CI-scale environments.
  expect_lt(r_ll$elapsed + r_lc$elapsed, 20)
})
