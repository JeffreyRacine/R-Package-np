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
  pilot <- get(".npreghat_complete")(bws = bw, txdat = x, exdat = x,
    y = y, output = "apply")
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
  # An explicit payload is authoritative, not recomputed from the bandwidth.
  args$fit.mean.train <- rep(mean(y), n)
  set.seed(812); first <- do.call(compute,args); seed <- .Random.seed
  set.seed(812); second <- do.call(compute,args)
  expect_identical(first[c("boot.err","bxp","boot.all.err")],
                   second[c("boot.err","bxp","boot.all.err")])
  expect_identical(seed,.Random.seed)
  args$fit.mean.train <- numeric(2L)
  expect_error(do.call(compute,args),"fit.mean.train payload is invalid")
}

test_that("wild regression missing pilots use explicit evaluation geometry", {
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
