test_that("conditional GNN preserves structural deleted-support admission", {
  old <- options(np.messages=FALSE, np.tree=FALSE)
  on.exit(options(old), add=TRUE)
  set.seed(2005)
  x <- rt(60,2); y <- .5*x+rt(60,3)
  for (tree in c(FALSE,TRUE)) for (degree in 2:3) {
    options(np.tree=tree)
    for (xx in list(x,x+10,3*x)) {
      X <- data.frame(x=xx); Y <- data.frame(y=y)
      bw <- npcdensbw(xdat=X, ydat=Y, bws=c(8,degree+1),
        bwtype="generalized_nn", bwmethod="cv.ls", regtype="lp",
        degree=degree, cxkertype="uniform", bandwidth.compute=FALSE)
      expect_identical(np:::.npcdensbw_eval_only(X,Y,bw,
        invalid.penalty="dbmax")$objective, -.Machine$double.xmax)
    }
  }
})

test_that("conditional GNN coefficient rejection returns an invalid trial and recovers", {
  old <- options(np.messages=FALSE, np.tree=FALSE)
  on.exit(options(old), add=TRUE)
  # Fixed input doubles from the independent 80-digit degree-9 diagnostic.
  X <- data.frame(x=c(-0.97664596000686288,
    -0.43895895127207041,
    -0.43034884426742792,
    -0.33993896842002869,
    -0.10370887443423271,
    -0.030859079211950302,
    0.060696707572788,
    0.097869373857975006,
    0.13481667498126626,
    0.15460366429761052,
    0.28428055951371789,
    0.31269944133237004,
    0.42187898280099034,
    0.48529501864686608,
    0.49828091962262988,
    0.52924963971599936,
    0.56254933495074511,
    0.5904632699675858,
    0.62544923089444637,
    0.66115597495809197,
    0.77588330907747149,
    0.89903982635587454,
    0.90697889775037766,
    0.98832042934373021))
  Y <- data.frame(y=c(-0.78128134739543742,
    -1.8117775577250161,
    -0.51426183918320545,
    -0.069188444154038242,
    -0.046298207433846966,
    -0.081255023005929433,
    -0.038633913380860613,
    0.79079310335597341,
    0.11373804769627469,
    0.14886092447722751,
    0.20678347055277324,
    0.59074023604058057,
    0.70308283891077938,
    1.1038119360494143,
    0.47845367472470157,
    1.3272642680891682,
    0.60417702829195252,
    1.2812540384856481,
    0.82652935750751011,
    1.9438440261005678,
    1.6689775446752155,
    1.5212322785244152,
    0.96193624295715918,
    0.61314715841749134))
  for(kernel in c("uniform", "epanechnikov")) for(tree in c(FALSE, TRUE)) {
    options(np.tree=tree)
    b <- npcdensbw(xdat=X, ydat=Y, bws=c(8,20), bwtype="generalized_nn",
      bwmethod="cv.ls", regtype="lp", degree=9L, bernstein.basis=FALSE,
      cxkertype=kernel, bandwidth.compute=FALSE)
    expect_identical(np:::.npcdensbw_eval_only(X,Y,b,
      invalid.penalty="dbmax")$objective, -.Machine$double.xmax)
    expect_identical(np:::.npcdensbw_eval_only(X,Y,b,
      invalid.penalty="baseline")$objective, -1e7)
  }
  # Search continues through invalid trials; an all-invalid final point is rejected.
  expect_error(npcdensbw(xdat=X, ydat=Y, bws=c(8,20),
    bwtype="generalized_nn", bwmethod="cv.ls", regtype="lp", degree=9L,
    bernstein.basis=FALSE, cxkertype="uniform", bwsolver="mads", nmulti=1L,
    nomad.opts=list(MAX_BB_EVAL=12L)),
    "did not return a raw-valid solution")
  expect_error(npcdensbw(xdat=X, ydat=Y, bws=c(8,20),
    bwtype="generalized_nn", bwmethod="cv.ls", regtype="lp", degree=9L,
    bernstein.basis=FALSE, cxkertype="epanechnikov", nmulti=1L),
    "kernel support or numerical accuracy", fixed=TRUE)
  # A stable degree-9 design keeps the existing public option available.
  X <- data.frame(x=sort(cos(pi*(0:59)/59)))
  Y <- data.frame(y=sin(2*X$x)+cos(3*X$x)/3)
  for(kernel in c("uniform", "epanechnikov")) {
    b <- npcdensbw(xdat=X, ydat=Y, bws=c(20,59), bwtype="generalized_nn",
      bwmethod="cv.ls", regtype="lp", degree=9L, bernstein.basis=FALSE,
      cxkertype=kernel, bandwidth.compute=FALSE)
    value <- np:::.npcdensbw_eval_only(X,Y,b, invalid.penalty="dbmax")$objective
    expect_true(is.finite(value) && abs(value) < .Machine$double.xmax)
  }
})

test_that("conditional GNN H60 retains its independently integrated objective", {
  old <- options(np.messages=FALSE, np.tree=FALSE)
  on.exit(options(old), add=TRUE)
  set.seed(2005);x <- rt(60,2);y <- .5*x+rt(60,3)
  X <- data.frame(x);Y <- data.frame(y)
  b <- npcdensbw(xdat=X, ydat=Y, bws=c(8,4), bwtype="generalized_nn",
    bwmethod="cv.ls", regtype="lp", degree=2L, bernstein.basis=FALSE,
    cxkertype="uniform", bandwidth.compute=FALSE)
  value <- np:::.npcdensbw_eval_only(X,Y,b)$objective
  reference <- -249132.23910049172
  expect_lt(abs(value-reference)/(1+abs(reference)),1e-10)
})
