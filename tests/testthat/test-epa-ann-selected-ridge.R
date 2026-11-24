test_that("Epanechnikov ANN hats report the ridge selected at each coordinate", {
  v <- readRDS(test_path("fixtures","r112-epanechnikov-ann.rds"))
  x <- data.frame(x=v$x)
  bw <- npregbw(xdat=x, ydat=v$y, bws=53, bandwidth.compute=FALSE,
    regtype="lp", degree=3, basis="glp", bernstein.basis=TRUE,
    bwtype="adaptive_nn", ckertype="epanechnikov", ckerorder=2)
  g <- data.frame(x=c(v$eval[2L], 1.1496080229039194))
  for (order in 0:1) {
    H <- npreghat(bw,txdat=x,exdat=g,s=order)
    ridge <- attr(H,"ridge.used",exact=TRUE)
    expect_equal(ridge,c(0,0.0003811790352210174),tolerance=1e-12)
    upper <- npreghat(bw,txdat=x,exdat=g[2L,,drop=FALSE],s=order)
    expect_identical(attr(upper,"ridge.used",exact=TRUE),ridge[2L])
  }
})
