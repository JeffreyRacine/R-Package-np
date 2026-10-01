r23_scale_spread <- function(z) min(sd(z), IQR(z)/(2*qnorm(.75)), mad(z))

test_that("conditional scale reference exercises the IQR branch", {
  # Sparse tails make sd large; asymmetric central spacing makes the
  # normal-consistent MAD exceed IQR/QFAC. No RNG or fragile near-tie.
  z <- rep(c(-100,-2,-1,-.5,0,3,3,4,100),each=5L)
  pieces <- c(sd(z),IQR(z)/(2*qnorm(.75)),mad(z))
  expect_identical(which.min(pieces),2L)
  expect_equal(r23_scale_spread(z),pieces[2L],tolerance=0)
  expect_gt(abs(r23_scale_spread(z)-IQR(z)/1.349),1e-6)
  pkg <- getNamespaceName(environment(npcdensbw))
  expect_equal(as.numeric(getFromNamespace("EssDee",pkg)(z)),
               r23_scale_spread(z),tolerance=2e-12)
})

test_that("conditional fixed scaling owns each response spread", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(9131)
  n <- 45L
  x <- data.frame(a = rnorm(n, sd = .2), b = runif(n, 0, 5))
  yy <- data.frame(y = 30*x$a + 2*x$b + rnorm(n, sd = 3),
                   z = -6*x$a + 3*x$b + rnorm(n, sd = .6))
  pkg <- getNamespaceName(environment(npcdensbw))
  for (family in c("density-ml", "density-ls", "distribution")) for (ny in 1:2) {
    cdf <- family == "distribution"
    make <- if (cdf) npcdistbw else npcdensbw
    fit <- if (cdf) npcdist else npcdens
    evaluate <- getFromNamespace(if (cdf) ".npcdistbw_eval_only" else ".npcdensbw_eval_only", pkg)
    y <- yy[seq_len(ny)]
    sf <- c(c(1, .7)[seq_len(ny)], .9, 1.2)
    args <- list(xdat = x, ydat = y, regtype = "lp", degree = c(1L, 1L),
                 bwmethod = if (family == "density-ml") "cv.ml" else "cv.ls",
                 bandwidth.compute = FALSE)
    scaled <- do.call(make, c(args, list(bws = sf, bwscaling = TRUE)))
    h <- sf*c(vapply(y, r23_scale_spread, 0.0), vapply(x, r23_scale_spread, 0.0))*
      n^(-1/(4+ncol(x)+ncol(y)))
    expect_equal(c(scaled$bandwidth$y, scaled$bandwidth$x), unname(h), tolerance = 2e-12)
    physical <- do.call(make, c(args, list(bws = unname(h), bwscaling = FALSE)))
    a <- evaluate(xdat = x, ydat = y, bws = scaled)$objective
    b <- evaluate(xdat = x, ydat = y, bws = physical)$objective
    expect_true(is.finite(a) && is.finite(b))
    expect_lte(abs(a-b), 2e-10*max(1,abs(a),abs(b)))
    args$xdat <- x[2:1]
    swapped <- do.call(make, c(args, list(bws = sf[c(seq_len(ny),ny+2L,ny+1L)],
                                         bwscaling = TRUE)))
    c <- evaluate(xdat = x[2:1], ydat = y, bws = swapped)$objective
    expect_lte(abs(a-c), 2e-10*max(1,abs(a),abs(c)))
    fs <- fitted(fit(scaled,txdat=x,tydat=y,se=FALSE))
    fp <- fitted(fit(physical,txdat=x,tydat=y,se=FALSE))
    expect_equal(fs,fp,tolerance=2e-10)
  }
})
