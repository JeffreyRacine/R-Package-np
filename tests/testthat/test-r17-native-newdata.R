test_that("native bandwidth newdata uses the same sample as explicit evaluation", {
  set.seed(17119)
  x <- data.frame(x = runif(32), z = runif(32))
  y <- data.frame(y = x$x + x$z + rnorm(32, sd = .2))
  e <- x[1:4, , drop = FALSE]
  ey <- y[1:4, , drop = FALSE]
  b <- npregbw(xdat = x, ydat = y$y, bws = c(.4, .4), bandwidth.compute = FALSE)
  expect_identical(fitted(npreg(b, newdata = e)), fitted(npreg(b, exdat = e)))
  expect_identical(fitted(predict(b, newdata = e)), fitted(npreg(b, exdat = e)))
  expect_length(fitted(npreg(b, newdata = e)), 4L)
  expect_identical(fitted(npreg(b, newdata = x, exdat = e)), fitted(npreg(b, exdat = e)))
  calls <- 0L
  a <- npreg(b, newdata = { calls <- calls + 1L; e })
  expect_identical(calls, 1L)
  expect_length(fitted(a), 4L)
  for (family in c("npudens", "npudist")) {
    b <- get(paste0(family, "bw"))(dat = x, bws = c(.4, .4), bandwidth.compute = FALSE)
    f <- get(family)
    expect_identical(fitted(f(b, newdata = e)), fitted(f(b, edat = e)))
    expect_identical(fitted(f(b, newdata = x, edat = e)), fitted(f(b, edat = e)))
  }
  for (family in c("npcdens", "npcdist")) {
    b <- get(paste0(family, "bw"))(xdat = x, ydat = y,
      bws = c(.4, .4, .4), bandwidth.compute = FALSE)
    f <- get(family)
    expect_identical(fitted(f(b, newdata = cbind(e, ey))), fitted(f(b, exdat = e, eydat = ey)))
    expect_identical(fitted(f(b, newdata = cbind(ey, e[, 2:1]))), fitted(f(b, exdat = e, eydat = ey)))
    expect_error(f(b, newdata = e), "must include columns")
    expect_identical(fitted(f(b, newdata = x, exdat = e, eydat = ey)), fitted(f(b, exdat = e, eydat = ey)))
  }
})

test_that("native newdata roles cover semiparametric and promoted estimators", {
  set.seed(18119)
  x <- data.frame(x = runif(40))
  z <- data.frame(z = runif(40))
  y <- x$x + sin(z$z) + rnorm(40, sd = .2)
  e <- x[1:4, , drop = FALSE]
  ez <- z[1:4, , drop = FALSE]
  b <- npindexbw(xdat = cbind(x, z), ydat = y, bws = c(.4, .5, 1),
                bandwidth.compute = FALSE)
  expect_identical(fitted(npindex(b, newdata = cbind(e, ez), se = FALSE)),
                   fitted(npindex(b, exdat = cbind(e, ez), se = FALSE)))
  for (family in c("npplreg", "npscoef")) {
    b <- get(paste0(family, "bw"))(xdat = x, ydat = y, zdat = z,
      bws = if (family == "npplreg") matrix(.4, 2L, 1L) else .4,
      bandwidth.compute = FALSE)
    f <- get(family)
    expect_identical(fitted(f(b, newdata = cbind(ez, e))),
                     fitted(f(b, exdat = e, ezdat = ez)))
    expect_length(fitted(f(b, newdata = cbind(e, ez))), 4L)
    expect_error(f(b, newdata = e), "must include columns")
    expect_identical(fitted(f(b, newdata = cbind(x, z), exdat = e, ezdat = ez)),
                     fitted(f(b, exdat = e, ezdat = ez)))
  }
  b <- npcdistbw(xdat = x, ydat = data.frame(y), bws = c(.4, .4),
                bandwidth.compute = FALSE)
  calls <- 0L
  a <- npqreg(b, newdata = { calls <- calls + 1L; e }, tau = c(.25, .5))
  expect_identical(calls, 1L)
  expect_identical(fitted(a), fitted(npqreg(b, exdat = e, tau = c(.25, .5))))
  expect_equal(dim(fitted(a)), c(4L, 2L))
  b <- npcdensbw(xdat = x, ydat = data.frame(y = factor(rep(0:1, 20))),
                bws = c(.4, .2), bandwidth.compute = FALSE)
  a <- npconmode(b, newdata = e, probabilities = TRUE)
  ref <- npconmode(b, exdat = e, probabilities = TRUE)
  expect_identical(fitted(a), fitted(ref))
  expect_length(fitted(a), 4L)
  expect_identical(predict(a, newdata = e, type = "prob"),
                   predict(ref, exdat = e, type = "prob"))
})
