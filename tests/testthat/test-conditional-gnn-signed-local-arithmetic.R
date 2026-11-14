# CF-210: frozen saved-data references from 256-bit original-design rows
# and independent direct response integration; no package route is the oracle.
r79_eval <- function(f, degree, order, kx, ky, pkg) {
  x <- data.frame(x=f$x); y <- data.frame(y=f$y)
  b <- npcdensbw(xdat=x,ydat=y,bws=c(ky,kx),bwtype="generalized_nn",
    bwmethod="cv.ls",regtype="lp",degree=degree,bernstein.basis=TRUE,
    cxkertype="epanechnikov",cxkerorder=order,bandwidth.compute=FALSE)
  snapshot <- serialize(b,NULL)
  value <- getFromNamespace(".npcdensbw_eval_only",pkg)(x,y,b,
    invalid.penalty="dbmax")$objective
  expect_identical(serialize(b,NULL),snapshot)
  value
}

test_that("signed GNN CVLS preserves the corrected row in its integral", {

  old <- options(np.messages=FALSE,np.tree=FALSE,np.extendednn=FALSE)
  on.exit(options(old),add=TRUE)
  f <- readRDS(test_path("fixtures","cf210_signed_gnn.rds"))
  reference <- c(0.0269694986456600,-0.0427000790059809,-0.140702022266686)
  for(j in seq_along(reference)) {
    got <- r79_eval(f$cluster,2L,c(4L,6L,8L)[j],12L,8L,"np")
    expect_lt(abs(got-reference[j]),1e-9*(1+abs(reference[j])))
  }
  # A high-degree control exposed rounded binomial factors after the row fix.
  for(o in c(4L,8L)) {
    ref <- if(o==4L)-3.72253885378886 else -2.15236612430578
    got <- r79_eval(f$ordinary,12L,o,35L,12L,"np")
    expect_lt(abs(got-ref),1e-9*(1+abs(ref)))
  }
  # Positive-kernel arithmetic and a pre-existing sparse refusal are retained.
  expect_equal(r79_eval(f$ordinary,2L,2L,35L,12L,"np"),
    0.759364178074412,tolerance=1e-12)
  expect_identical(r79_eval(f$ordinary,3L,4L,2L,12L,"np"),
    -.Machine$double.xmax)
})
