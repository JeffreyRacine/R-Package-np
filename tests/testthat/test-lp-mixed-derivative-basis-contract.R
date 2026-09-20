test_that("LP derivative products annihilate every differentiated constant factor", {
  set.seed(920311)
  x <- data.frame(a=runif(80,-1,1),b=runif(80,2,4))
  e <- data.frame(a=c(-.7,0,.4,.8),b=c(2.2,3,3.4,3.8))
  for (basis in c("glp","additive","tensor")) for (bernstein in c(FALSE,TRUE)) {
    degree <- c(2L,2L)
    terms <- npBuildLpTerms(degree,basis)
    W <- W.lp(x,degree=degree,basis=basis,bernstein.basis=bernstein)
    raw <- W.lp(x,degree=degree,basis=basis,bernstein.basis=FALSE)
    # Exact change of basis; all monomials, not just a selected fitted response.
    coef <- qr.solve(W,raw)
    for (order in list(c(0L,0L),c(1L,0L),c(2L,0L),c(1L,1L),c(2L,1L))) {
      D <- W.lp(x,exdat=e,degree=degree,gradient.vec=order,
                basis=basis,bernstein.basis=bernstein)
      oracle <- vapply(seq_len(nrow(terms)),function(j) {
        if (any(terms[j,] < order)) return(rep(0,nrow(e)))
        prod(factorial(terms[j,])/factorial(terms[j,]-order)) *
          e$a^(terms[j,1]-order[1]) * e$b^(terms[j,2]-order[2])
      },numeric(nrow(e)))
      expect_equal(unname(D %*% coef),oracle,tolerance=1e-10)
      expect_equal(unname(W.lp(x,exdat=e[1,,drop=FALSE],degree=degree,
        gradient.vec=order,basis=basis,bernstein.basis=bernstein)),
        unname(D[1,,drop=FALSE]),tolerance=0)
    }
    expect_equal(W.lp(x,exdat=e,degree=degree,gradient.vec=c(0,0),
      basis=basis,bernstein.basis=bernstein),
      W.lp(x,exdat=e,degree=degree,basis=basis,bernstein.basis=bernstein),
      tolerance=0)
  }
  expect_equal(W.lp(x,degree=c(0,0),gradient.vec=c(1,0)),matrix(0,80,1))
  expect_equal(unname(W.lp(x,degree=c(0,2),gradient.vec=c(1,0))),matrix(0,80,3))
})

test_that("mixed regression hats reproduce a polynomial cross derivative", {
  old <- options(np.messages=FALSE)
  on.exit(options(old),add=TRUE)
  set.seed(920311)
  x <- data.frame(a=runif(80,-1,1),b=runif(80,2,4))
  e <- data.frame(a=c(-.7,0,.4,.8),b=c(2.2,3,3.4,3.8))
  y <- with(x,2+3*a+2*b+4*a^2+5*a*b)
  for (bernstein in c(FALSE,TRUE)) {
    bw <- npregbw(xdat=x,ydat=y,bws=c(.7,.7),bandwidth.compute=FALSE,
                  regtype="lp",degree=c(2,2),bernstein.basis=bernstein)
    H <- npreghat(bw,txdat=x,exdat=e,s=c(1,1),output="matrix")
    expect_equal(drop(H %*% y),rep(5,nrow(e)),tolerance=1e-10)
    expect_equal(as.double(npreghat(bw,txdat=x,exdat=e,y=y,s=c(1,1),output="apply")),
                 rep(5,nrow(e)),tolerance=1e-10)
  }
})
