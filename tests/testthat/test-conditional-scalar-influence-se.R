test_that("conditional beta response allows a zero explanatory derivative total", {
  x <- data.frame(x = seq(.05, .95, length.out = 12L))
  y <- data.frame(y = c(.12, .65, .31, .82, .24, .57, .43, .91, .18, .73, .36, .52))
  for (cdf in c(FALSE, TRUE)) {
    constructor <- if (cdf) npcdistbw else npcdensbw
    estimator <- if (cdf) npcdist else npcdens
    b <- constructor(xdat = x, ydat = y, bws = c(.15, 2),
      bandwidth.compute = FALSE, regtype = "lc", cxkertype = "uniform",
      cykertype = "beta", cykerbound = "fixed", cykerlb = 0, cykerub = 1)
    fit <- estimator(bws = b, txdat = x, tydat = y, exdat = data.frame(x = .5),
      eydat = data.frame(y = .43), gradients = TRUE)
    z <- as.vector(npksum(bws = .15, txdat = y, exdat = data.frame(y = .43),
      ckertype = "beta", ckerbound = "fixed", ckerlb = 0, ckerub = 1,
      operator = if (cdf) "integral" else "normal", return.kernel.weights = TRUE)$kw)
    expect_equal(as.vector(fitted(fit)), mean(z), tolerance = 2e-12)
    expect_identical(as.vector(gradients(fit)), 0)
    expect_identical(as.vector(fit$congerr), 0)
    # Invalid base normalization is checked in the serial counterpart.
    # Do not exercise the separately deferred native-unwind pool cleanup here.
  }
})
