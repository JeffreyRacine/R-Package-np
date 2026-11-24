test_that("beta matrices retain the chosen ridge at each evaluation coordinate", {
  v <- readRDS(test_path("fixtures", "r112-cubic-beta-ann.rds"))
  x <- data.frame(x=v$x)
  bw <- npregbw(xdat=x, ydat=v$y, bws=24, bandwidth.compute=FALSE,
    regtype="lp", degree=3, basis="glp", bernstein.basis=TRUE,
    bwtype="adaptive_nn", ckertype="beta", ckerorder=2,
    ckerbound="fixed", ckerlb=-1, ckerub=1)
  ex <- data.frame(x=v$eval[c(2L,4L)])
  for (s in 0:1) {
    H <- npreghat(bw, txdat=x, exdat=ex, s=s, output="matrix")
    ridge <- attr(H,"ridge.used",exact=TRUE)
    expect_equal(ridge[1],0)
    expect_equal(ridge[2],0.017775181958286507,tolerance=1e-12)
    single <- npreghat(bw,txdat=x,exdat=ex[2,,drop=FALSE],s=s,output="matrix")
    expect_equal(attr(single,"ridge.used",exact=TRUE),ridge[2])
    expect_equal(as.double(single),as.double(H[2,]),tolerance=1e-12)
  }
})
