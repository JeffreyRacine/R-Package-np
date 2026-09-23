test_that("hat owners use one resolved polynomial specification", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(861)
  x <- data.frame(x = runif(45, .1, .9), z = runif(45, .1, .9))
  y <- sin(4 * x$x) + x$z^2 + x$x * x$z
  e <- data.frame(x = c(.22, .4, .65, .8), z = c(.7, .6, .35, .2))
  requests <- list(
    list(degree = c(1L, 1L), basis = "glp", bernstein.basis = FALSE),
    list(degree = c(1L, 1L), basis = "tensor", bernstein.basis = FALSE),
    list(degree = c(1L, 1L), basis = "additive", bernstein.basis = TRUE))
  for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
    for (kernel in c("gaussian", "epanechnikov", "beta")) {
      args <- list(xdat = x, ydat = y, bws = rep(if (bwtype == "fixed") .25 else 18, 2),
                   bwtype = bwtype, ckertype = kernel, regtype = "lp",
                   degree = c(2L, 2L), bandwidth.compute = FALSE)
      if (kernel == "beta")
        args <- c(args, list(ckerbound = "fixed", ckerlb = 0, ckerub = 1))
      b <- do.call(npregbw, args)
      original <- serialize(b, NULL)
      for (request in requests) {
        requested.args <- args
        requested.args[names(request)] <- request
        control <- do.call(npregbw, requested.args)
        for (s in 0:1) {
          call <- c(list(bws = b, txdat = x, exdat = e, s = s), request)
          H <- do.call(npreghat, call)
          ref <- npreghat(control, txdat = x, exdat = e, s = s)
          a <- do.call(npreghat, c(call, list(output = "apply", y = y)))
          many <- do.call(npreghat, c(call, list(output = "apply", y = cbind(y, 2 * y))))
          constrained <- do.call(npreghat, c(call, list(output = "constraint", y = y)))
          expect_equal(as.vector(H), as.vector(ref), tolerance = 1e-9)
          expect_equal(as.vector(a), as.vector(ref %*% y), tolerance = 1e-9)
          expect_equal(as.vector(many[, 1L]), as.vector(a), tolerance = 1e-9)
          expect_equal(as.vector(many[, 2L]), 2 * as.vector(a), tolerance = 1e-9)
          expect_equal(as.vector(constrained), as.vector(t(ref) * y), tolerance = 1e-9)
          expect_equal(as.vector(predict(H, exdat = e, y = y, output = "apply")),
                       as.vector(a), tolerance = 1e-9)
          expect_identical(attr(H, "bws")$degree.engine, control$degree.engine)
          expect_identical(attr(H, "bws")$basis.engine, control$basis.engine)
          expect_identical(attr(H, "bws")$bernstein.basis.engine,
                           control$bernstein.basis.engine)
        }
      }
      expect_identical(serialize(b, NULL), original)
    }
  }
})

test_that("higher-degree generalized hat fallback also retains overrides", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = seq(.1, .9, length.out = 23))
  y <- sin(4 * x$x)
  e <- data.frame(x = c(.2, .4, .7))
  b <- npregbw(xdat = x, ydat = y, bws = 9, bwtype = "generalized_nn",
               regtype = "lp", degree = 3L, bandwidth.compute = FALSE)
  control <- npregbw(xdat = x, ydat = y, bws = 9, bwtype = "generalized_nn",
                     regtype = "lp", degree = 2L, bandwidth.compute = FALSE)
  H <- npreghat(b, txdat = x, exdat = e, degree = 2L, s = 1L)
  ref <- npreghat(control, txdat = x, exdat = e, s = 1L)
  expect_equal(as.vector(H), as.vector(ref), tolerance = 1e-9)
  expect_equal(as.vector(npreghat(b, txdat = x, exdat = e, degree = 2L,
                                  s = 1L, output = "apply", y = y)),
               as.vector(ref %*% y), tolerance = 1e-9)
})
