test_that("factor, copula and hat native predictions use the same name rule", {
  if (exists("spawn_mpi_slaves", mode = "function")) {
    spawn_mpi_slaves()
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }
  set.seed(1941)
  n <- 36L
  x <- data.frame(a = runif(n), u = factor(rep(c("a", "b"), n/2)),
                  o = ordered(rep(0:2, n/3)))
  y <- x$a + as.numeric(x$u) + as.numeric(x$o)/4
  b <- npregbw(xdat = x, ydat = y, bws = c(.5, .2, .2),
               bandwidth.compute = FALSE, regtype = "ll")
  e <- x[c(4, 8, 10), ]
  a <- npreg(b, se = FALSE)
  expect_equal(predict(a, newdata = e[, 3:1]), predict(a, exdat = e))
  expect_equal(predict(a, newdata = e[1, 3:1]), predict(a, exdat = e[1, ]))
  expect_error(predict(a, newdata = e[c("a", "u")]), "must include columns")
  du <- e; names(du) <- c("a", "a", "o")
  expect_error(predict(a, newdata = du), "ambiguous")
  H <- npreghat(b, txdat = x)
  expect_equal(as.numeric(predict(H, newdata = e[, 3:1], output = "apply", y = y)),
               as.numeric(predict(H, exdat = e, output = "apply", y = y)))
  expect_equal(as.numeric(predict(H, newdata = data.frame(wrong = 0),
                                 exdat = e, output = "apply", y = y)),
               as.numeric(predict(H, exdat = e, output = "apply", y = y)))

  mode.bw <- npcdensbw(xdat = x, ydat = data.frame(g = x$u),
                       bws = c(.5, .2, .2, .2), bandwidth.compute = FALSE)
  mode <- npconmode(mode.bw, probabilities = TRUE)
  expect_equal(predict(mode, newdata = e[, 3:1]), predict(mode, exdat = e))
  expect_equal(fitted(npconmode(mode.bw, newdata = e[, 3:1])),
               fitted(npconmode(mode.bw, exdat = e)))

  cx <- data.frame(a = runif(n), b = runif(n))
  cb <- npudistbw(dat = cx, bws = c(.3, .3), bandwidth.compute = FALSE)
  co <- npcopula(cb, data = cx, u = c(.3, .7), se = FALSE,
                 n.quasi.inv = 20, er.quasi.inv = .2)
  u <- data.frame(u1 = c(.3, .4), u2 = c(.6, .7))
  expect_equal(predict(co, newdata = u[, 2:1]), predict(co, u = u))
  expect_error(predict(co, newdata = data.frame(a = .3, b = .7)),
               "probability coordinates")
  du <- u; names(du) <- c("u1", "u1")
  expect_error(predict(co, newdata = du), "ambiguous")
})

test_that("the native newdata role splitter preserves named and positional contracts", {
  d <- data.frame(x = 1:3, z = 4:6)
  split <- function(nd) .np_native_newdata_parts(nd, list(exdat = c("x", "z")), "test")$exdat
  expect_identical(split(d[, 2:1]), d)
  expect_identical(split(cbind(d[, 2:1], extra = 7)), d)
  expect_equal(split(as.matrix(d[, 2:1])), d)
  expect_equal(unname(as.matrix(split(unname(as.matrix(d))))), unname(as.matrix(d)))
  expect_error(split(d["x"]), "must include columns")
  dup <- d; names(dup) <- c("x", "x")
  expect_error(split(dup), "ambiguous")
  expect_identical(dim(split(d[1L, 2:1])), c(1L, 2L))
})

test_that("native fit and prediction entry points match newdata by name", {
  skip_on_cran()
  if (exists("spawn_mpi_slaves", mode = "function")) {
    spawn_mpi_slaves()
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }
  set.seed(91841)
  n <- 36L
  x <- data.frame(a = runif(n, -1, 1), b = runif(n, -1, 1))
  z <- data.frame(z = runif(n, -1, 1))
  y <- 2*x$a - x$b + sin(z$z) + rnorm(n, sd = .2)
  yd <- data.frame(y = y)
  e <- x[c(2, 8, 20), ]; ez <- z[c(2, 8, 20), , drop = FALSE]
  ey <- yd[c(2, 8, 20), , drop = FALSE]
  for (fam in c("npreg", "npudens", "npudist", "npcdens", "npcdist",
                "npqreg", "npplreg", "npscoef", "npindex", "nplsqreg")) {
    ctor <- get(paste0(if (fam == "npqreg") "npcdist" else fam, "bw"))
    fit <- get(fam)
    u <- fam %in% c("npudens", "npudist")
    cond <- fam %in% c("npcdens", "npcdist", "npqreg")
    semi <- fam %in% c("npplreg", "npscoef")
    ba <- if (u) list(dat = x, bws = c(.45, .45)) else if (cond)
      list(xdat = x, ydat = yd, bws = c(.45, .45, .45)) else if (semi)
      list(xdat = x, ydat = y, zdat = z,
           bws = if (fam == "npplreg") matrix(.45, 3, 1) else .45) else if (fam == "npindex")
      list(xdat = x, ydat = y, bws = c(.45, 1, .5)) else
      list(xdat = x, ydat = y, bws = c(.45, .45))
    ba$bandwidth.compute <- FALSE
    if (fam %in% c("npreg", "nplsqreg")) ba$regtype <- "ll"
    b <- do.call(ctor, ba)
    a <- fit(bws = b, se = FALSE)
    nd <- if (semi) cbind(e, ez) else if (cond && fam != "npqreg") cbind(e, ey) else e
    native <- if (u) list(edat = e) else if (semi) list(exdat = e, ezdat = ez) else
      if (cond && fam != "npqreg") list(exdat = e, eydat = ey) else list(exdat = e)
    oracle <- fitted(do.call(fit, c(list(bws = b, se = FALSE), native)))
    shuffled <- nd[, rev(seq_len(ncol(nd))), drop = FALSE]
    expect_equal(fitted(fit(bws = b, newdata = shuffled, se = FALSE)), oracle,
                 tolerance = 1e-10, info = fam)
    expect_equal(predict(a, newdata = shuffled), oracle, tolerance = 1e-10, info = fam)
    expect_equal(predict(a, newdata = cbind(shuffled, extra = 0)), oracle,
                 tolerance = 1e-10, info = fam)
    expect_equal(do.call(predict, c(list(object = a, newdata = data.frame(wrong = 0)), native)),
                 oracle, tolerance = 1e-10, info = fam)
  }
})
