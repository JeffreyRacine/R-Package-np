test_that("bounded check and KS losses use the directed deleted fit", {
  ns <- asNamespace(getNamespaceName(environment(npregbw)))
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(4040)
  n <- 60L
  x <- runif(n, -1, 1)
  y <- sin(2*x) + rnorm(n, sd = .3)
  yb <- as.numeric(y > median(y))
  X <- data.frame(x = x)
  check <- get(".nplsqreg_call_fixed_degree_core", ns)
  evaluate <- get(".npregbw_eval_only", ns)
  for (degree in 0:2) {
    # Removing a positive query-only boundary mass leaves this WLS fit intact.
    # This reference does not call a package fitting or objective helper.
    predictions <- function(response) vapply(seq_len(n), function(i) {
      design <- outer(x[-i]-x[i], 0:degree, `^`)
      sw <- sqrt(dnorm((x[-i]-x[i])/.3)/.3)
      q <- qr(design*sw)
      stopifnot(q$rank == degree+1L)
      qr.coef(q, response[-i]*sw)[1L]
    }, 0)
    fitted <- predictions(y)
    probability <- pmax(sqrt(.Machine$double.eps),
      pmin(1-sqrt(.Machine$double.eps), predictions(yb)))
    expected.ks <- -mean(yb*log(probability)+(1-yb)*log1p(-probability))
    for (tree in c(FALSE, TRUE)) {
      options(np.tree = tree)
      bw <- npregbw(xdat=X, ydat=y, bws=.3, regtype="lp", degree=degree,
        ckerbound="fixed", ckerlb=-1, ckerub=1, bandwidth.compute=FALSE)
      loss <- function() check(X,y,rep(1,n),.5,bw,.5,c(.01,.99),
                              list(invalid.penalty="dbmax"),FALSE)$objective
      got <- if (getNamespaceName(ns) == "npRmpi")
        get(".npRmpi_with_local_regression", ns)(loss()) else loss()
      expect_equal(got, mean(abs(y-fitted))/2, tolerance=1e-10)
      expect_equal(evaluate(X,yb,bw,objective="ks")$objective,
                   expected.ks, tolerance=1e-10)
      expect_equal(evaluate(X,y,bw)$objective,
                   mean((y-fitted)^2), tolerance=1e-10)
    }
  }
})
