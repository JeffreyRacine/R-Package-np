rank_support_fixture <- function() {
  set.seed(4711)
  x <- data.frame(x = runif(40, .05, .95), z = runif(40, .05, .95))
  list(x = x, y = sin(4*x$x) + x$z^2 + rnorm(40, sd = .2),
       ex = x[c(3, 9, 14, 22, 31), ])
}

test_that("structurally sparse LL hats and fits share the canonical ridge", {
  pkg <- environmentName(environment(npreg))
  f <- rank_support_fixture()
  bw <- suppressWarnings(npregbw(xdat = f$x, ydat = f$y, bws = c(11, 11),
    bandwidth.compute = FALSE, bwtype = "generalized_nn", ckertype = "uniform",
    regtype = "ll"))
  kw <- getFromNamespace(".np_kernel_weights_direct", pkg)(
    bws = bw, txdat = f$x, exdat = f$ex, bandwidth.divide = TRUE)
  basis <- getFromNamespace("W.lp", pkg)
  z <- basis(f$x, degree = c(1L, 1L), basis = "glp", bernstein.basis = FALSE)
  d <- basis(f$x, exdat = f$ex, degree = c(1L, 1L), basis = "glp",
             bernstein.basis = FALSE, gradient.vec = c(1L, 0L))
  expect_identical(sum(kw[, 1L] != 0), 2L)
  a <- crossprod(z, z * kw[, 1L])
  delta <- max(abs(diag(a))) / nrow(z)
  u <- solve(a + diag(delta, ncol(z)), d[1L, ])
  u[1L] <- u[1L] * (1 + delta / a[1L, 1L])
  oracle <- kw[, 1L] * drop(z %*% u)
  old <- getOption("matprod")
  on.exit(options(matprod = old), add = TRUE)
  for (mode in c("default", "blas", "internal")) {
    options(matprod = mode)
    h <- npreghat(bw, txdat = f$x, exdat = f$ex, s = c(1L, 0L),
                  output = "matrix")
    expect_equal(as.double(h[1L, ]), oracle, tolerance = 2e-12, info = mode)
    for (y in list(f$y, cos(f$x$x), rep(1, 40))) {
      g <- gradients(npreg(bw, txdat = f$x, tydat = y, exdat = f$ex,
                            gradients = TRUE))[, 1L]
      expect_equal(drop(h %*% y), g, tolerance = 2e-10, info = mode)
    }
  }
})

test_that("donor certificates count rows rather than bootstrap multiplicity", {
  pkg <- environmentName(environment(npreg))
  count <- getFromNamespace(".np_inid_lp_rank_bounds", pkg)
  project <- getFromNamespace(".np_inid_lp_batch_project", pkg)
  w <- c(0, 1, -1, 0, 2)
  counts <- cbind(c(0, 8, 0, 0, 0), c(0, 1, 1, 0, 6), rep(1, 5))
  expect_identical(count(w, 3L, counts), c(1L, 3L, 3L))
  for (k in 0:4)
    expect_identical(count(c(rep(1, k), rep(0, 5-k)), 3L), as.integer(min(3, k)))
  z <- cbind(1, seq(-1, 1, length.out = 5), seq(-1, 1, length.out = 5)^2)
  y <- sin(seq_len(5))
  cw <- w * counts[, 1L]
  a <- crossprod(z, z*cw)
  rhs <- drop(crossprod(z, y*cw))
  packed <- unlist(lapply(1:3, function(j) a[j, j:3]))
  result <- project(matrix(packed, 1), matrix(rhs, 1), c(1, .2, .04),
                    represented.mass = 8, diagnostics = TRUE, rank.bounds = 1L)
  expect_identical(result$ridge_steps, 1L)
  delta <- max(abs(diag(a))) / 8
  rhs[1L] <- rhs[1L] * (1 + delta/a[1L, 1L])
  expected <- sum(c(1, .2, .04) * solve(a + diag(delta, 3), rhs))
  expect_equal(result$values[1L, 1L], expected, tolerance = 2e-13)
  expect_error(project(matrix(packed, 1), matrix(rhs, 1), c(1, 0, 0),
                       represented.mass = 8, rank.bounds = -1L), "rank_bounds")
})
test_that("higher-degree sparse rows use the same response and adjoint policy", {
  pkg <- environmentName(environment(npreg))
  x <- data.frame(x = c(-.1, 0, .1, seq(1, 3, length.out=9)),
                  z = c(.1, -.1, 0, 2.3, 1.1, 2.7, 1.4, 3, 1.8, 2.1, 1.3, 2.6))
  y <- sin(x$x) + cos(x$z)
  ex <- data.frame(x = 0, z = 0)
  for (basis.name in c("glp", "additive", "tensor"))
  for (bernstein in c(FALSE, TRUE)) {
    bw <- npregbw(xdat=x, ydat=y, bws=c(.3,.3), regtype="lp",
      degree=c(2L,2L), degree.select="manual", basis=basis.name,
      bernstein.basis=bernstein,
      ckertype="uniform", bandwidth.compute=FALSE)
    kw <- getFromNamespace(".np_kernel_weights_direct",pkg)(
      bws=bw, txdat=x, exdat=ex, bandwidth.divide=TRUE)[,1L]
    basis <- getFromNamespace("W.lp",pkg)
    z <- basis(x,degree=c(2L,2L),basis=basis.name,bernstein.basis=bernstein)
    expect_lt(sum(kw!=0),ncol(z))
    a <- crossprod(z,z*kw)
    delta <- max(abs(diag(a)))/nrow(z)
    for (s in list(c(0L,0L),c(1L,0L),c(0L,1L))) {
      d <- basis(x,exdat=ex,degree=c(2L,2L),basis=basis.name,
        bernstein.basis=bernstein,gradient.vec=if(sum(s)) s else NULL)
      u <- solve(a+diag(delta,ncol(z)),d[1L,])
      u[1L] <- u[1L]*(1+delta/a[1L,1L])
      oracle <- kw*drop(z%*%u)
      h <- npreghat(bw,txdat=x,exdat=ex,s=s,output="matrix")
      expect_equal(as.double(h),oracle,tolerance=2e-11)
      fit <- npreg(bw,txdat=x,tydat=y,exdat=ex,gradients=TRUE)
      answer <- if(sum(s)==0L) fitted(fit) else gradients(fit)[,which(s!=0L)]
      expect_equal(drop(h%*%y),as.double(answer),tolerance=2e-11)
    }
  }
})
