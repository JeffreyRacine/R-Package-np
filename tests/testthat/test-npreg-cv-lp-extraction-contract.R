library(np)

make_npreg_cv_lp_data <- function(n = 120L, seed = 20260305L) {
  set.seed(seed)
  x1 <- runif(n)
  x2 <- runif(n)
  z1 <- factor(rbinom(n, 1L, 0.5))
  z2 <- ordered(rbinom(n, 1L, 0.5))
  y <- cos(2 * pi * x1) + 0.5 * sin(2 * pi * x2) + as.numeric(z1) + rnorm(n, sd = 0.18)
  data.frame(y = y, x1 = x1, x2 = x2, z1 = z1, z2 = z2)
}

make_npreg_cv_lp_cont_data <- function(n = 120L, seed = 20260305L) {
  set.seed(seed)
  x <- runif(n)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.12)
  data.frame(y = y, x = x)
}

run_npreg_cv_lp_case <- function(dat, formula, regtype, tree, degree = NULL, basis = NULL) {
  old_opts <- options(np.messages = FALSE, np.tree = tree)
  on.exit(options(old_opts), add = TRUE)

  args <- list(
    formula = formula,
    data = dat,
    regtype = regtype,
    bwmethod = "cv.ls",
    nmulti = 1
  )
  if (!is.null(degree))
    args$degree <- degree
  if (!is.null(basis))
    args$basis <- basis

  do.call(npregbw, args)
}

test_that("npregbw extracted LP CV path preserves ll == lp degree-1 under tree on/off", {
  dat <- make_npreg_cv_lp_cont_data()

  for (tree in c(FALSE, TRUE)) {
    set.seed(404)
    bw.ll <- run_npreg_cv_lp_case(dat, formula = y ~ x, regtype = "ll", tree = tree)
    set.seed(404)
    bw.lp <- run_npreg_cv_lp_case(dat, formula = y ~ x, regtype = "lp", tree = tree, degree = 1L, basis = "glp")

    expect_equal(as.numeric(bw.ll$fval), as.numeric(bw.lp$fval), tolerance = 1e-12)
    expect_equal(as.numeric(bw.ll$bw), as.numeric(bw.lp$bw), tolerance = 1e-10)
    expect_equal(as.integer(bw.ll$num.feval), as.integer(bw.lp$num.feval))
  }
})

test_that("npregbw extracted LP CV path is deterministic and tree-invariant for degree-1 GLP", {
  dat <- make_npreg_cv_lp_data(seed = 20260306L)

  set.seed(505)
  bw.f0.a <- run_npreg_cv_lp_case(dat, formula = y ~ x1 + x2 + z1 + z2, regtype = "lp", tree = FALSE, degree = c(1L, 1L), basis = "glp")
  set.seed(505)
  bw.f0.b <- run_npreg_cv_lp_case(dat, formula = y ~ x1 + x2 + z1 + z2, regtype = "lp", tree = FALSE, degree = c(1L, 1L), basis = "glp")
  set.seed(505)
  bw.t1 <- run_npreg_cv_lp_case(dat, formula = y ~ x1 + x2 + z1 + z2, regtype = "lp", tree = TRUE, degree = c(1L, 1L), basis = "glp")
  fit.f0 <- npreg(bws = bw.f0.a, bandwidth.compute = FALSE)
  fit.t1 <- npreg(bws = bw.t1, bandwidth.compute = FALSE)

  expect_equal(as.numeric(bw.f0.a$fval), as.numeric(bw.f0.b$fval), tolerance = 1e-12)
  expect_equal(as.numeric(bw.f0.a$bw), as.numeric(bw.f0.b$bw), tolerance = 1e-12)
  expect_lt(abs(as.numeric(bw.f0.a$fval) - as.numeric(bw.t1$fval)), 2e-9)
  expect_equal(fitted(fit.f0), fitted(fit.t1), tolerance = 1e-4)
})
