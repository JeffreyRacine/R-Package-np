test_that("private derivative-error demand preserves existing request masks", {
  expect_identical(.np_regression_output_request(FALSE, FALSE), 0L)
  expect_identical(.np_regression_output_request(TRUE, FALSE), 1L)
  expect_identical(.np_regression_output_request(FALSE, TRUE), 2L)
  expect_identical(.np_regression_output_request(TRUE, TRUE), 3L)
  expect_identical(.np_regression_output_request(TRUE, TRUE, FALSE), 7L)
  expect_error(.np_regression_output_request(TRUE, TRUE, NA), "TRUE or FALSE")
  expect_error(.np_regression_output_request(FALSE, TRUE, TRUE), "require both")
  expect_error(.np_regression_output_request(TRUE, FALSE, TRUE), "require both")
})

test_that("mean uncertainty and derivative values do not require derivative errors", {
  if (exists("spawn_mpi_slaves", mode = "function") && !spawn_mpi_slaves())
    skip("Could not spawn MPI slaves")
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(928L)
  x <- data.frame(x = runif(96L, -1, 1))
  y <- sin(2*x$x) + rnorm(96L, sd = .1)
  for (rt in c("lc", "ll", "lp")) {
    bw.args <- list(xdat=x, ydat=y, bws=.45, bandwidth.compute=FALSE,
                    regtype=rt)
    if (rt == "lp")
      bw.args <- c(bw.args, list(degree=2L, bernstein.basis=TRUE))
    b <- do.call(npregbw, bw.args)
    full <- npreg(bws=b, txdat=x, tydat=y, se=TRUE, gradients=TRUE)
    lean <- npreg(bws=b, txdat=x, tydat=y, se=TRUE, gradients=TRUE,
                  .np.gradient.errors=FALSE)
    expect_identical(lean$mean, full$mean, info=rt)
    expect_identical(lean$merr, full$merr, info=rt)
    expect_identical(lean$grad, full$grad, info=rt)
    expect_null(lean$gerr)
    expect_true(all(is.finite(full$gerr)))
    expect_error(npreg(bws=b, txdat=x, tydat=y, se=FALSE, gradients=TRUE,
                      .np.gradient.errors=TRUE), "require both")
  }
})

test_that("default index covariance remains available without public gradients", {
  if (exists("spawn_mpi_slaves", mode = "function") && !spawn_mpi_slaves())
    skip("Could not spawn MPI slaves")
  old <- options(np.messages=FALSE, np.tree=FALSE)
  on.exit(options(old), add=TRUE)
  set.seed(9281L)
  d <- data.frame(x=rnorm(128L), z=rnorm(128L))
  for (method in c("kleinspady", "ichimura")) {
    d$y <- if (method == "kleinspady")
      as.numeric(d$x + .6*d$z + rnorm(128L) > 0)
    else sin(d$x + .6*d$z) + rnorm(128L, sd=.2)
    b <- npindexbw(y~x+z, data=d, method=method, regtype="ll",
                   bws=c(1,.6,.6), bandwidth.compute=FALSE)
    lean <- npindex(bws=b, data=d)
    full <- npindex(bws=b, data=d, gradients=TRUE)
    expect_true(lean$se)
    expect_identical(coef(lean), coef(full))
    expect_identical(vcov(lean), vcov(full))
    expect_identical(lean$mean, full$mean)
    expect_identical(lean$merr, full$merr)
    expect_true(is.finite((coef(lean)/sqrt(diag(vcov(lean))))[2L]))
    expect_null(lean$gerr)
    expect_true(all(is.finite(full$gerr)))
  }
})

test_that("beta mean uncertainty does not require gradient uncertainty", {
  if (exists("spawn_mpi_slaves", mode = "function") && !spawn_mpi_slaves())
    skip("Could not spawn MPI slaves")
  old <- options(np.messages=FALSE, np.tree=FALSE)
  on.exit(options(old), add=TRUE)
  set.seed(9286L)
  x <- data.frame(x=runif(96L,.1,.9), u=factor(rep(letters[1:3],32L)),
                  o=ordered(rep(1:3,32L)))
  y <- sin(2*x$x) + rnorm(96L,sd=.1)
  for (mixed in c(FALSE,TRUE)) for (bt in c("fixed","generalized_nn","adaptive_nn")) {
    xx <- if (mixed) x else x["x"]
    h <- if (bt=="fixed") .14 else 64
    b <- npregbw(xdat=xx,ydat=y,bws=c(h,if(mixed)c(.25,.25)),
      bwtype=bt,ckertype="beta",ckerbound="fixed",ckerlb=0,ckerub=1,
      regtype="lc",bandwidth.compute=FALSE)
    full <- npreg(bws=b,txdat=xx,tydat=y,se=TRUE,gradients=TRUE)
    lean <- npreg(bws=b,txdat=xx,tydat=y,se=TRUE,gradients=TRUE,
                  .np.gradient.errors=FALSE)
    expect_identical(lean$mean,full$mean)
    expect_identical(lean$merr,full$merr)
    expect_identical(lean$grad,full$grad)
    expect_null(lean$gerr)
    expect_true(all(is.finite(full$gerr)))
    ex <- xx[1:7,,drop=FALSE]
    ex$x <- seq(.2,.8,length.out=7L)
    full.external <- npreg(bws=b,txdat=xx,tydat=y,exdat=ex,
                           se=TRUE,gradients=TRUE)
    lean.external <- npreg(bws=b,txdat=xx,tydat=y,exdat=ex,
                           se=TRUE,gradients=TRUE,.np.gradient.errors=FALSE)
    expect_identical(lean.external$mean,full.external$mean)
    expect_identical(lean.external$merr,full.external$merr)
    expect_identical(lean.external$grad,full.external$grad)
    expect_null(lean.external$gerr)
  }
})
