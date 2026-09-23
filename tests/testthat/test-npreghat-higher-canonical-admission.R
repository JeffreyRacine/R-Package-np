test_that("translated higher and mixed hats retain canonical polynomial admission", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(9283)
  x <- data.frame(a = runif(61, -1, 1), b = runif(61, -1, 1))
  e <- data.frame(a = c(-.4, .2, .5), b = c(.3, -.2, .7))
  y <- 1 + 2*x$a + 3*x$b + 4*x$a*x$b + 5*x$a^2
  for (shift in c(0, 60)) {
    xx <- x + shift; ee <- e + shift
    b <- npregbw(xdat = xx, ydat = y, bws = c(.7, .8), regtype = "lp",
      degree = c(2, 2), basis = "glp", bernstein.basis = FALSE,
      bandwidth.compute = FALSE)
    for (s in list(c(2, 0), c(1, 1))) {
      H <- npreghat(b, txdat = xx, exdat = ee, s = s)
      truth <- rep(if (s[1] == 2) 10 else 4, 3)
      # Absolute error allows raw polynomial Gram rounding after translation;
      # the pre-repair unwanted ridge produces errors of order one instead.
      expect_lt(max(abs(drop(H %*% y) - truth)), 2e-5)
      expect_identical(attr(H, "ridge.used"), rep(0, 3))
      a <- npreghat(b, txdat = xx, exdat = ee, s = s,
                    y = cbind(y, 2*y), output = "apply")
      expect_equal(unname(a), unname(H %*% cbind(y, 2*y)), tolerance = 1e-12)
    }
  }
})

test_that("higher hats share canonical prepared solves across basis and geometry", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  pkg <- environmentName(environment(npregbw))
  ns <- asNamespace(pkg)
  weights <- get(".np_kernel_weights_direct", ns)
  design <- get("W.lp", ns)
  set.seed(8252)
  x <- data.frame(a = runif(37, -.8, .8), b = runif(37, -.8, .8))
  e <- x[c(3, 8, 15), ] + .01
  y <- sin(x$a) + x$b^2
  for (type in c("fixed", "generalized_nn", "adaptive_nn"))
    for (basis in c("glp", "tensor")) for (bernstein in c(FALSE, TRUE))
      for (tree in c(FALSE, TRUE)) {
        options(np.tree = tree)
        b <- npregbw(xdat = x, ydat = y,
          bws = if (type == "fixed") c(.7, .8) else c(27, 29),
          bwtype = type, regtype = "lp", degree = c(2, 2), basis = basis,
          bernstein.basis = bernstein, bandwidth.compute = FALSE)
        for (s in list(c(2L, 0L), c(1L, 1L))) {
          H <- npreghat(b, txdat = x, exdat = e, s = s)
          K <- weights(b, x, e)
          W <- design(x, degree = c(2L, 2L), basis = basis,
                      bernstein.basis = bernstein)
          E <- design(x, exdat = e, degree = c(2L, 2L), gradient.vec = s,
                      basis = basis, bernstein.basis = bernstein)
          oracle <- t(vapply(seq_len(nrow(e)), function(i) {
            v <- .Call("C_np_lp_adjoint_prepared", crossprod(W, W*K[, i]),
              as.double(E[i, ]), as.integer(nrow(x)),
              as.integer(min(ncol(W), sum(K[, i] != 0))), PACKAGE = pkg)
            K[, i] * drop(W %*% v)
          }, numeric(nrow(x))))
          expect_equal(as.vector(H), as.vector(oracle), tolerance = 1e-10)
          expect_true(all(is.finite(attr(H, "ridge.used"))))
          constraint <- npreghat(b, txdat = x, exdat = e, s = s,
                                y = y, output = "constraint")
          expect_equal(as.vector(constraint), as.vector(t(H)*y), tolerance = 1e-12)
        }
      }
})

test_that("bounded and associated-kernel higher hats retain polynomial reproduction", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = seq(.1, .9, length.out = 41))
  e <- data.frame(x = c(.22, .43, .7))
  y <- x$x^2
  for (type in c("fixed", "generalized_nn", "adaptive_nn"))
    for (kernel in c("gaussian", "epanechnikov", "beta"))
      for (bernstein in c(FALSE, TRUE)) {
        b <- npregbw(xdat = x, ydat = y, bws = if (type == "fixed") .6 else 29,
          bwtype = type, regtype = "lp", degree = 2, ckertype = kernel,
          ckerbound = "fixed", ckerlb = 0, ckerub = 1,
          bernstein.basis = bernstein, bandwidth.compute = FALSE)
        H <- npreghat(b, txdat = x, exdat = e, s = 2)
        expect_equal(drop(H %*% y), rep(2, nrow(e)), tolerance = 1e-10)
        expect_identical(attr(H, "ridge.used"), rep(0, nrow(e)))
        expect_equal(as.vector(predict(H, exdat = e, y = y, output = "apply")),
                     rep(2, nrow(e)), tolerance = 1e-10)
      }
})

test_that("explicit higher-derivative ridge keeps its requested absolute operator", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  ns <- asNamespace(environmentName(environment(npregbw)))
  x <- data.frame(x = c(-1, -.7, -.2, .1, .4, .8, 1))
  e <- data.frame(x = c(-.3, .6))
  y <- sin(x$x)
  b <- npregbw(xdat = x, ydat = y, bws = .6, regtype = "lp", degree = 2,
    bernstein.basis = FALSE, bandwidth.compute = FALSE)
  K <- get(".np_kernel_weights_direct", ns)(b, x, e)
  W <- cbind(1, x$x, x$x^2)
  for (rho in c(.1, 2)) {
    H <- npreghat(b, txdat = x, exdat = e, s = 2, ridge = rho)
    oracle <- t(vapply(seq_len(nrow(e)), function(i) {
      A <- crossprod(W, W * K[, i]) + diag(rho, ncol(W))
      K[, i] * drop(W %*% solve(t(A), c(0, 0, 2)))
    }, numeric(nrow(x))))
    expect_equal(as.vector(H), as.vector(oracle), tolerance = 1e-12)
    expect_identical(attr(H, "ridge.used"), rep(rho, nrow(e)))
  }
})

test_that("default higher hats preserve undefined external rows without fabricating zeros", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = seq(-1, 1, length.out = 17))
  e <- data.frame(x = c(.1, 9))
  b <- npregbw(xdat = x, ydat = x$x^2, bws = .7, regtype = "lp", degree = 2,
    ckertype = "epanechnikov", bandwidth.compute = FALSE)
  expect_warning(H <- npreghat(b, txdat = x, exdat = e, s = 2),
                 "all computed kernel weights are zero")
  expect_true(all(is.finite(H[1, ])))
  expect_true(all(is.na(H[2, ])))
  expect_true(is.na(attr(H, "ridge.used")[2]))
  expect_warning(a <- npreghat(b, txdat = x, exdat = e, s = 2,
                               y = x$x^2, output = "apply"),
                 "all computed kernel weights are zero")
  expect_equal(a[1], 2, tolerance = 1e-10)
  expect_true(is.na(a[2]))
  expect_warning(c <- npreghat(b, txdat = x, exdat = e, s = 2,
                               y = x$x^2, output = "constraint"),
                 "all computed kernel weights are zero")
  expect_true(all(is.na(c[, 2])))
  expect_error(npreghat(b, txdat = x, exdat = e, s = 2,
                       .np.require.finite = TRUE), "LP solve failed")
  H.ridge <- npreghat(b, txdat = x, exdat = e, s = 2, ridge = .1)
  expect_identical(unname(H.ridge[2, ]), rep(0, nrow(x)))
  expect_identical(attr(H.ridge, "ridge.used"), rep(.1, nrow(e)))
})
