test_that("conditional leading SE uses separate X and Y kernel moments", {
  n <- 41L
  x <- data.frame(x1 = seq(-.4, .4, length.out = n),
                  x2 = seq(-.2, .2, length.out = n))
  y <- data.frame(y = .4 * sin(seq_len(n)))
  ex <- data.frame(x1 = c(-.03, .02), x2 = c(.01, -.02))
  ey <- data.frame(y = c(.02, -.01))
  moments <- c(gaussian = 1/(2*sqrt(pi)), uniform = .5)
  for (cdf in c(FALSE, TRUE)) for (kx in names(moments))
    for (ky in names(moments)) {
      bwfun <- if(cdf) npcdistbw else npcdensbw
      fitfun <- if(cdf) npcdist else npcdens
      bw <- bwfun(xdat=x, ydat=y, bws=rep(.8,3),
                  bandwidth.compute=FALSE, regtype="lc",
                  cxkertype=kx, cykertype=ky)
      fit <- fitfun(bws=bw,txdat=x,tydat=y,exdat=ex,eydat=ey)
      wx <- vapply(seq_len(nrow(ex)), function(q) {
        products <- lapply(seq_len(ncol(x)), function(j) {
          u <- (x[[j]]-ex[[j]][q])/.8
          if(kx=="gaussian") dnorm(u) else .5*(abs(u)<1)
        })
        sum(Reduce(`*`,products))
      },numeric(1))
      value <- as.double(fitted(fit))
      variance <- if(cdf) value*(1-value)*moments[[kx]]^2/wx else
        value*moments[[kx]]^2*moments[[ky]]/(.8*wx)
      expect_equal(as.double(se(fit)),sqrt(variance),tolerance=2e-12)
    }
})
