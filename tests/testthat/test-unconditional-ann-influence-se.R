local({
  n <- 48L
  i <- seq_len(n)
  x <- data.frame(a = sin(i*.714)+i/300, b = cos(i*.417)+sin(i*.71)/3,
                  o = ordered(rep(1:3, length.out=n)))
  at <- x[c(7L, 23L, 39L), ]; at$a <- at$a+.013

  # Dense diagnostic only: independent rank geometry and explicit interval
  # matrix. Canonical kernel weights distinguish geometry from kernel tests.
  oracle <- function(bw, density) {
    kbw <- kbandwidth(bw)
    # Independent explicit density-family kernel mapping, checked below
    # against the actual public point, not just the SE producer's adapter.
    if (identical(bw$okertype, "liracine")) kbw$okertype <- "nliracine"
    kw <- npksum(bws=kbw, txdat=x, exdat=at,
      operator=if(density) "normal" else "integral",
      permutation.operator=if(density) "derivative" else "normal",
      bandwidth.divide=TRUE, return.kernel.weights=TRUE,
      return.derivative.kernel.weights=TRUE,
      .np.internal.bandwidth.divide.weights=TRUE)
    ff <- if (density) npudens else npudist
    expect_equal(colMeans(kw$kw), fitted(ff(bws=bw,tdat=x,edat=at,se=FALSE)),
                 tolerance=2e-12)
    k <- c(11L, 15L)
    phi <- kw$kw
    for (j in seq_len(2L)) {
      z <- x[[j]]
      distance <- abs(outer(z,z,"-")); diag(distance) <- Inf
      ranks <- apply(distance,1L,sort)
      m <- min(floor((n-1)^(2/3)), floor(min(k[j],n-1-k[j])/2))
      r <- ranks[k[j], ]; b <- 2*m/(n-1)/(ranks[k[j]+m, ]-ranks[k[j]-m, ])
      inside <- abs(outer(z,z,"-")) <= r
      for (q in seq_len(nrow(at))) {
        w <- ((at[[j]][q]-z)*kw$p.kw[,q,j] +
              if(density) kw$kw[,q] else 0)/(r*b)
        phi[,q] <- phi[,q] + as.vector(crossprod(inside,w))/n
      }
    }
    sqrt(apply(phi,2L,var)/n)
  }

  run <- function() {
    old <- options(np.messages=FALSE, np.tree=getOption("np.tree", "auto"),
                   np.extendednn=TRUE)
    on.exit(options(old), add=TRUE)
    for (density in c(TRUE,FALSE)) for (kernel in c("gaussian","epanechnikov")) {
      test_that(paste("unconditional ANN joint influence",density,kernel), {
        bf <- if(density) npudensbw else npudistbw
        ff <- if(density) npudens else npudist
        bw <- bf(dat=x,bws=c(11,15,.3),bwtype="adaptive_nn",
                 ckertype=kernel,bandwidth.compute=FALSE)
        for (tree in list(FALSE,TRUE,"auto")) {
          options(np.tree=tree)
          off <- ff(bws=bw,tdat=x,edat=at,se=FALSE)
          on <- ff(bws=bw,tdat=x,edat=at,se=TRUE)
          expect_identical(fitted(off),fitted(on))
          expect_length(off$derr,0L)
          expect_equal(se(on),oracle(bw,density),tolerance=2e-11)
        }
      })
    }
    test_that("unqualified ANN uncertainty keeps points and says why", {
      bw <- npudensbw(dat=x,bws=c(2*n,2*n,.3),bwtype="adaptive_nn",bandwidth.compute=FALSE)
      expect_warning(on <- npudens(bws=bw,tdat=x,edat=at,se=TRUE),
                     "no regular two-sided interior spacing pilot")
      expect_true(all(is.na(se(on))))
      off <- npudens(bws=bw,tdat=x,edat=at,se=FALSE)
      expect_identical(fitted(on),fitted(off))
    })
    for (density in c(TRUE, FALSE)) for (ordered in c("wangvanryzin", "liracine", "racineliyan")) {
      test_that(paste("ANN covariance matches actual ordered kernel",density,ordered), {
        bf <- if(density) npudensbw else npudistbw
        ff <- if(density) npudens else npudist
        bw <- bf(dat=x,bws=c(11,15,.3),bwtype="adaptive_nn",okertype=ordered,
                 bandwidth.compute=FALSE)
        on <- ff(bws=bw,tdat=x,edat=at,se=TRUE)
        expect_equal(se(on),oracle(bw,density),tolerance=2e-11)
      })
    }
  }
  if (exists(".npRmpi_with_local_regression", mode="function"))
    .npRmpi_with_local_regression(run())
  else run()
})
