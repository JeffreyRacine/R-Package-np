b1_wild_case <- function(type, regtype, degree = NULL, basis = "glp",
                         bernstein = FALSE, tree = FALSE) {
  old <- options(np.messages = FALSE, np.tree = tree)
  on.exit(options(old), add = TRUE)
  set.seed(42)
  n <- 60L
  x <- data.frame(x = runif(n), v = rnorm(n),
    u = factor(sample(letters[1:3], n, TRUE)),
    o = ordered(sample(1:3, n, TRUE)))
  y <- sin(3*x$x) + .2*x$v + .3*(x$u == "b") +
    .1*as.integer(x$o) + rnorm(n, sd = .2)
  bw <- npregbw(xdat = x, ydat = y, bws = c(
    if (type == "fixed") c(.25,.6) else c(15,18), .3,.4),
    bandwidth.compute = FALSE, bwtype = type, regtype = regtype,
    degree = degree, basis = basis, bernstein.basis = bernstein)
  get <- function(name) getFromNamespace(name, "npRmpi")
  # Independent public fitted-row target, without an external self-query.
  pilot <- fitted(npreg(bws = bw, txdat = x, tydat = y))
  if (identical(type, "generalized_nn")) {
    external <- fitted(npreg(bws = bw, txdat = x, tydat = y, exdat = x))
    expect_gt(max(abs(pilot - external)), 1e-4)
  }
  compute <- get("compute.bootstrap.errors.rbandwidth")
  ex <- x[c(7,20,40),,drop = FALSE]
  # A categorical panel is a common frame with exactly one row per level.
  cat.ex <- ex[rep(1L,3L),,drop = FALSE]
  cat.ex$u <- factor(letters[1:3], levels = levels(x$u))
  for (mode in c("level","continuous","categorical")) {
    args <- list(xdat = x, ydat = y, exdat = if (mode=="categorical") cat.ex else ex,
      gradients = mode!="level", gradient.order = 1L,
      slice.index = if (mode=="categorical") 3L else 1L,
      plot.errors.boot.method = "wild", plot.errors.boot.nonfixed = "exact",
      plot.errors.boot.wild = "rademacher", plot.errors.boot.blocklen = 1L,
      plot.errors.boot.num = 7L, plot.errors.center = "estimate",
      plot.errors.type = "pmzsd", plot.errors.alpha = .05, bws = bw)
    set.seed(811)
    actual <- do.call(compute,args)
    actual.seed <- .Random.seed
    set.seed(811)
    expected <- do.call(compute,c(args,list(fit.mean.train=as.vector(pilot))))
    fields <- c("boot.err","bxp","boot.all.err")
    expect_equal(actual[fields],expected[fields],tolerance=5e-10)
    expect_identical(actual.seed,.Random.seed)
  }
  args$fit.mean.train <- numeric(2L)
  expect_error(do.call(compute,args),"fit.mean.train payload is invalid")
}

b1_supplied_pilot_case <- function(type,
    compute = getFromNamespace("compute.bootstrap.errors.rbandwidth", "npRmpi")) {
  withr::local_preserve_seed()
  withr::local_options(list(np.messages = FALSE, np.tree = FALSE))
  set.seed(13218)
  n <- 24L
  x <- data.frame(x = runif(n), u = factor(rep(letters[1:3], 8L)))
  y <- sin(2*x$x) + .2*(x$u == "b") + rnorm(n, sd = .2)
  width <- if (type == "fixed") .28 else 12L
  bw <- npregbw(xdat = x, ydat = y, bws = c(width, .3),
    bandwidth.compute = FALSE, bwtype = type, regtype = "ll",
    ckertype = "gaussian", ckerorder = 2L, ukertype = "aitchisonaitken")
  ex <- x[c(3L, 11L, 20L), , drop = FALSE]
  # Independent external-query LL rows; no package hat/bootstrap oracle.
  H <- t(vapply(seq_len(nrow(ex)), function(i) {
    dx <- x$x - ex$x[i]
    h <- if (type == "fixed") width else sort(abs(dx))[width]
    w <- dnorm(dx/h) * ifelse(x$u == ex$u[i], .7, .15)
    design <- cbind(1, dx)
    gram <- crossprod(design, w * design)
    expect_gt(rcond(gram), 1e-8)
    w * as.vector(design %*% solve(gram, c(1, 0)))
  }, numeric(n)))
  pilots <- list(.3*sin(seq_len(n)), -.7 + .4*cos(seq_len(n)/2))
  errors <- vector("list", length(pilots))
  for (i in seq_along(pilots)) {
    set.seed(1129)
    signs <- matrix(ifelse(runif(n*7L) <= .5, -1, 1), n, 7L)
    expected.seed <- .Random.seed
    draws <- t(H %*% (pilots[[i]] + (y - pilots[[i]]) * signs))
    expected <- qnorm(.975) * apply(draws, 2L, sd)
    set.seed(1129)
    actual <- compute(xdat = x, ydat = y, exdat = ex,
      fit.mean.train = pilots[[i]], gradients = FALSE, gradient.order = 1L,
      slice.index = 1L, plot.errors.boot.method = "wild",
      plot.errors.boot.nonfixed = "exact", plot.errors.boot.wild = "rademacher",
      plot.errors.boot.blocklen = 1L, plot.errors.boot.num = 7L,
      plot.errors.center = "estimate", plot.errors.type = "pmzsd",
      plot.errors.alpha = .05, bws = bw)
    expect_equal(actual$boot.err[, 1L], expected, tolerance = 5e-10)
    expect_equal(actual$boot.err[, 2L], expected, tolerance = 5e-10)
    expect_identical(.Random.seed, expected.seed)
    errors[[i]] <- actual$boot.err[, 1:2, drop = FALSE]
  }
  expect_false(identical(errors[[1L]], errors[[2L]]))
}

test_that("wild regression automatic pilots use fitted training geometry", {
  for (tree in c(FALSE,TRUE))
    for (type in c("fixed","generalized_nn","adaptive_nn"))
      for (regtype in c("lc","ll"))
        b1_wild_case(type,regtype,tree=tree)
})

test_that("wild regression pilot contract covers higher-degree bases", {
  for (basis in c("glp","additive","tensor"))
    for (bernstein in c(FALSE,TRUE))
      b1_wild_case("generalized_nn","lp",degree=c(2L,2L),basis=basis,
                   bernstein=bernstein)
})

test_that("supplied wild pilots determine the independent residual transform", {
  for (type in c("fixed", "generalized_nn"))
    b1_supplied_pilot_case(type)
})

test_that("exact wild derivative defaults match each public response refit", {
  withr::local_preserve_seed()
  withr::local_options(list(np.messages = FALSE, np.tree = FALSE))
  set.seed(711)
  n <- 32L
  x <- data.frame(x = runif(n), u = factor(rep(letters[1:2], n/2)))
  y <- sin(3*x$x) + .2*(x$u == "b") + rnorm(n, sd = .2)
  ex <- x[c(4L, 15L, 26L), , drop = FALSE]
  helper <- getFromNamespace(".np_wild_boot_from_regression_exact", "npRmpi")
  for (type in c("fixed", "generalized_nn", "adaptive_nn")) {
    bw <- npregbw(xdat = x, ydat = y,
      bws = c(if (type == "fixed") .3 else 14L, .25),
      bandwidth.compute = FALSE, bwtype = type, regtype = "lc")
    pilot <- fitted(npreg(bws = bw, txdat = x, tydat = y))
    set.seed(123)
    signs <- matrix(ifelse(runif(n*7L) <= .5, -1, 1), n, 7L)
    expected.rng <- .Random.seed
    expected <- t(vapply(seq_len(7L), function(j) {
      response <- pilot + (y-pilot)*signs[,j]
      gradients(npreg(bws = bw, txdat = x, tydat = response,
                      exdat = ex, gradients = TRUE))[,1L]
    }, numeric(nrow(ex))))
    set.seed(123)
    actual <- helper(xdat = x, exdat = ex, bws = bw, ydat = y,
                     B = 7L, wild = "rademacher", gradients = TRUE,
                     slice.index = 1L)
    expect_equal(actual$t, expected, tolerance = 5e-10)
    expect_identical(.Random.seed, expected.rng)
  }
})

test_that("public fit and bandwidth plots share an automatic training pilot", {
  withr::local_preserve_seed()
  withr::local_options(list(np.messages = FALSE, np.tree = FALSE))
  set.seed(719)
  x <- data.frame(x = runif(40), u = factor(rep(letters[1:2], 20)),
                  o = ordered(rep(1:2, each = 20)))
  y <- sin(3*x$x) + .2*(x$u == "b") + rnorm(40, sd = .2)
  fields <- c("mean", "grad", "merr", "gerr", "bxp")
  for (regtype in c("lc", "ll", "lp")) {
    bw <- npregbw(xdat = x, ydat = y, bws = c(16, .25, .3),
      bandwidth.compute = FALSE, bwtype = "generalized_nn",
      regtype = regtype, degree = if (regtype == "lp") 2L else NULL)
    model <- npreg(bws = bw, txdat = x, tydat = y)
    for (grad in c(FALSE, TRUE)) {
      set.seed(981)
      from.fit <- plot(model, output = "data", errors = "bootstrap",
        plot.errors.boot.method = "wild", B = 7L, neval = 5L,
        gradients = grad)
      fit.rng <- .Random.seed
      set.seed(981)
      from.bw <- plot(bw, xdat = x, ydat = y, output = "data",
        errors = "bootstrap", plot.errors.boot.method = "wild",
        B = 7L, neval = 5L, gradients = grad)
      expect_identical(.Random.seed, fit.rng)
      expect_equal(lapply(from.fit, function(z) z[fields]),
                   lapply(from.bw, function(z) z[fields]), tolerance = 5e-10)
    }
  }
})
