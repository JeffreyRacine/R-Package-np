.dc1_alllarge_fixture <- function(cdf, ridge) {
  n <- 24L
  tx <- data.frame(x1 = seq(-.6, .6, length.out = n))
  tx$x2 <- tx$x1
  if (!ridge) tx$x2[12L] <- .3
  ty <- data.frame(y = .5*sin(seq_len(n)*1.3) + seq_len(n)/60)
  ex <- data.frame(x1 = .1, x2 = .1)
  ey <- data.frame(y = .2)
  hy <- .25
  bw <- do.call(if (cdf) npcdistbw else npcdensbw,
    list(xdat = tx, ydat = ty, bws = c(hy, 1.3, 1.3),
         bandwidth.compute = FALSE, bwscaling = FALSE, bwtype = "fixed",
         regtype = "lp", degree = c(1L, 1L), basis = "glp",
         bernstein.basis = FALSE, cxkertype = "uniform", cykertype = "gaussian"))

  # Independent accepted all-large point map: no GENERAL intercept correction.
  # The CDF preserves its historical Gaussian-CDF constant.
  B <- cbind(1, tx$x1, tx$x2)
  A0 <- crossprod(B)
  z <- if (cdf) pnorm(sqrt(2)*.7071067810*(ey$y-ty$y)/hy) else
    dnorm((ey$y-ty$y)/hy)/hy
  r <- if (ridge) 1/n else 0
  C <- solve(A0 + diag(r, ncol(B)))
  d <- cbind(c(1, .1, .1), c(0, 1, 0), c(0, 0, 1))
  a <- t(C) %*% d
  sigma2 <- mean((z-mean(z))^2)
  list(args = list(bws = bw, txdat = tx, tydat = ty, exdat = ex, eydat = ey),
       point = drop(crossprod(d, C %*% crossprod(B, z))),
       error = sqrt(sigma2*diag(crossprod(a, A0 %*% a))),
       old.error = sqrt(sigma2*diag(crossprod(d, C %*% d))))
}

test_that("conditional all-large errors follow the accepted regularized map", {
  if (!spawn_mpi_slaves(1L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  for (cdf in c(FALSE, TRUE)) for (ridge in c(FALSE, TRUE)) {
    d <- .dc1_alllarge_fixture(cdf, ridge)
    fun <- if (cdf) npcdist else npcdens
    fit <- expect_warning(do.call(fun,
      c(d$args, list(gradients = TRUE, se = TRUE))), NA)
    off <- expect_warning(do.call(fun,
      c(d$args, list(gradients = TRUE, se = FALSE))), NA)
    expect_equal(c(as.numeric(fitted(fit)), as.numeric(gradients(fit))),
                 d$point, tolerance = 2e-11)
    expect_equal(c(as.numeric(se(fit)), as.numeric(gradients(fit, se = TRUE))),
                 d$error, tolerance = 2e-11)
    expect_identical(fitted(fit), fitted(off))
    expect_identical(gradients(fit), gradients(off))
    if (ridge) expect_gt(max(abs(d$error-d$old.error)), .1)
  }
})
