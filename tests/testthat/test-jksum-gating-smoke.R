library(npRmpi)

session_pool_active <- function(comm = 1L) {
  if (!isTRUE(getOption("npRmpi.mpi.initialized", FALSE)))
    return(FALSE)
  size <- tryCatch(as.integer(mpi.comm.size(comm)), error = function(e) NA_integer_)
  !is.na(size) && size >= 2L
}

ensure_session_slave_pool <- function() {
  if (!session_pool_active()) {
    suppressWarnings(npRmpi.init(nslaves = 1, quiet = TRUE))
  }
  invisible(TRUE)
}

with_session_slave_pool <- function(expr) {
  ensure_session_slave_pool()
  force(expr)
}

teardown({
  if (session_pool_active())
    try(npRmpi.quit(mode = "spawn", force = TRUE), silent = TRUE)
})

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

run_reg_cv_once <- function(dat, regtype, bwmethod, ...) {
  t_bw <- system.time(
    bw <- npregbw(
      y ~ x1 + x2 + z1 + z2,
      regtype = regtype,
      bwmethod = bwmethod,
      nmulti = 1,
      data = dat,
      ...
    )
  )

  list(
    fval = as.numeric(bw$fval),
    nfe = as.integer(bw$num.feval),
    nfe_fast = as.integer(bw$num.feval.fast),
    elapsed = as.numeric(t_bw[["elapsed"]])
  )
}

test_that("jksum regression CV parity is deterministic for mixed data", {
  skip_on_cran()
  old_opts <- options(
    np.messages = FALSE,
    np.largeh.rel.tol = 0.05,
    np.disc.upper.rel.tol = 0.05
  )
  on.exit(options(old_opts), add = TRUE)

  with_session_slave_pool({
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
})

test_that("jksum regression CV smoke performance remains bounded", {
  skip_on_cran()
  old_opts <- options(
    np.messages = FALSE,
    np.largeh.rel.tol = 0.05,
    np.disc.upper.rel.tol = 0.05
  )
  on.exit(options(old_opts), add = TRUE)

  with_session_slave_pool({
    dat <- make_jksum_mixed_data(n = 220L, seed = 101L)

    set.seed(321)
    r_ll <- run_reg_cv_once(dat, "ll", "cv.ls")
    set.seed(321)
    r_lc <- run_reg_cv_once(dat, "lc", "cv.ls")

    expect_true(r_ll$nfe > 0L)
    expect_true(r_lc$nfe > 0L)
    expect_true(is.finite(r_ll$fval))
    expect_true(is.finite(r_lc$fval))
    expect_lt(r_ll$elapsed + r_lc$elapsed, 20)
  })
})

test_that("large-h fast gateway is active again for lc/ll/lp under canonical DGPs", {
  skip_on_cran()
  old_opts <- options(
    np.messages = FALSE,
    np.largeh.rel.tol = 0.05,
    np.disc.upper.rel.tol = 0.05
  )
  on.exit(options(old_opts), add = TRUE)

  with_session_slave_pool({
    set.seed(42)
    n <- 200L
    x <- runif(n)
    y_lc <- rnorm(n, sd = 0.5 * sd(x))
    y_ll <- x + rnorm(n, sd = 0.5 * sd(x))
    dat_lc <- data.frame(y = y_lc, x = x, z1 = factor(0L), z2 = ordered(0L))
    dat_ll <- data.frame(y = y_ll, x = x, z1 = factor(0L), z2 = ordered(0L))

    set.seed(42)
    bw_ll <- npregbw(y ~ x, data = dat_ll, regtype = "ll", bwmethod = "cv.ls", nmulti = 1)
    set.seed(42)
    bw_lc <- npregbw(y ~ x, data = dat_lc, regtype = "lc", bwmethod = "cv.ls", nmulti = 1)
    set.seed(42)
    bw_lp <- npregbw(
      y ~ x,
      data = dat_ll,
      regtype = "lp",
      ckerorder = 4,
      bwmethod = "cv.ls",
      nmulti = 1
    )

    expect_true(is.finite(as.numeric(bw_ll$fval)))
    expect_true(is.finite(as.numeric(bw_lc$fval)))
    expect_true(is.finite(as.numeric(bw_lp$fval)))

    expect_gt(as.integer(bw_ll$num.feval.fast), 0L)
    expect_gt(as.integer(bw_lp$num.feval.fast), 0L)
    expect_gt(as.integer(bw_lc$num.feval.fast), 0L)
  })
})
