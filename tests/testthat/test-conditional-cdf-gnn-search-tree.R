test_that("GNN CDF search trees do not alter adjacent selection policies", {
  ns <- asNamespace("npRmpi")
  choose <- get(".npcdistbw_tree_code", ns)
  yes <- get("DO_TREE_YES", ns); no <- get("DO_TREE_NO", ns)
  old <- options(np.tree = FALSE, np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = seq(.1, .9, length.out = 14))
  y <- data.frame(y = sin(x$x))
  template <- npcdistbw(xdat = x, ydat = y, bws = c(7, 7),
    bwtype = "generalized_nn", cxkertype = "epanechnikov",
    cykertype = "epanechnikov", bandwidth.compute = FALSE)
  for (type in c("fixed", "generalized_nn", "adaptive_nn"))
    for (kernel in c("gaussian", "epanechnikov", "uniform", "beta"))
      for (mode in list(FALSE, TRUE, "auto")) for (context in c(FALSE, TRUE)) {
        b <- template; b$type <- type; b$cxkertype <- b$cykertype <- kernel
        options(np.tree = mode)
        compact <- kernel %in% c("epanechnikov", "uniform")
        expected <- isTRUE(mode) || (identical(mode, "auto") && compact)
        if (type == "generalized_nn") expected <- context && isTRUE(mode) && compact
        expect_identical(choose(b, 2L, 0L, cv.context = context), if (expected) yes else no)
      }
  options(np.tree = TRUE)
  for (role in c("cxkerorder", "cykerorder")) for (order in c(4L, 6L, 8L)) {
    b <- template; b[[role]] <- order
    expect_identical(choose(b, 2L, 0L, cv.context = TRUE), no)
  }
  for (role in c("cxkerbound", "cykerbound")) for (bound in c("range", "fixed")) {
    b <- template; b[[role]] <- bound
    expect_identical(choose(b, 2L, 0L, cv.context = TRUE), no)
  }
})

test_that("compact GNN CDF search trees preserve literal deleted-fold objectives", {
  old <- options(np.tree = TRUE, np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = c(.03,.03,.15,.21,.29,.36,.44,.51,.59,.67,.73,.81,.87,.93))
  y <- data.frame(y = c(.13,.13,.34,.22,.54,.26,.46,.65,.43,.73,.64,.82,.71,.91))
  n <- nrow(x); k <- 7L
  for (kernel in c("epanechnikov", "uniform")) for (degree in 0:2)
    for (empirical in c(FALSE, TRUE)) {
      K <- if (kernel == "uniform") function(z) .5*(abs(z)<1) else
        function(z) pmax(0, 3/(4*sqrt(5))*(1-z*z/5))
      P <- if (kernel == "uniform") function(z) pmin(1,pmax(0,.5+.5*z)) else
        function(z) {q <- pmax(-sqrt(5),pmin(sqrt(5),z)); .5+3/(4*sqrt(5))*(q-q^3/15)}
      grid <- if (empirical) y else y[c(14,7,1,4),,drop=FALSE]
      expected <- 0
      for (i in seq_len(n)) {
        keep <- setdiff(seq_len(n),i)
        hx <- sort(abs(x$x[keep]-x$x[i]))[k]
        z <- (x$x[keep]-x$x[i])/hx
        design <- outer(z,0:degree,"^"); w <- K(z)
        influence <- as.vector(design %*% solve(crossprod(design,w*design),
          c(1,rep(0,degree))))*w
        for (j in seq_len(nrow(grid))) {
          if (empirical && i==j) next
          donor <- if (empirical) setdiff(keep,j) else keep
          hy <- sort(abs(y$y[donor]-grid$y[j]))[k]
          prediction <- sum(influence*P((grid$y[j]-y$y[keep])/hy))
          expected <- expected+(as.integer(y$y[i]<=grid$y[j])-prediction)^2
        }
      }
      expected <- expected/(n*(nrow(grid)-as.integer(empirical)))
      b <- npcdistbw(xdat=x,ydat=y,bws=c(k,k),bandwidth.compute=FALSE,
        bwtype="generalized_nn",regtype="lp",degree=degree,
        cxkertype=kernel,cykertype=kernel)
      observed <- getFromNamespace(".npcdistbw_eval_only","npRmpi")(
        xdat=x,ydat=y,bws=b,gydat=if (empirical) NULL else grid,
        do.full.integral=empirical,invalid.penalty="dbmax")$objective
      expect_true(is.finite(observed) && abs(observed-expected)<2e-10,
        info=paste(kernel,degree,empirical,observed,expected))
    }
})
