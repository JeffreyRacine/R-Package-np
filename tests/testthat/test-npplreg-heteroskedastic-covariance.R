test_that("partial-linear covariance retains paired heteroskedastic scores", {
  ns <- asNamespace(getNamespaceName(environment(npplreg)))
  solve.linear <- get(".np_plreg_linear_solve",ns)
  g <- expand.grid(r=c(-2,-1,1,2),s=c(-2,-1,1,2),e=c(-1,1))
  X <- cbind(first=g$r,second=.3*g$r+g$s)
  u <- (g$r+.5*g$s)*g$e
  yhat <- seq_len(nrow(X))*.01
  y <- as.vector(X%*%c(2,-.7)+u)
  for(columns in list(1L,1:2,2:1)) {
    xx <- X[,columns,drop=FALSE]
    yy <- if(length(columns)==1L) 2*xx[,1L]+u else y
    fit <- solve.linear(yy,xx,yy+yhat,yhat,zdim=1L)
    inv <- solve(crossprod(xx))
    error <- yy+yhat-fit$train.fit
    meat <- Reduce(`+`,lapply(seq_len(nrow(xx)),function(i)
      error[i]^2*tcrossprod(xx[i,])))
    target <- nrow(xx)/(nrow(xx)-ncol(xx)-1L)*inv%*%meat%*%inv
    expect_equal(unname(fit$vcov),unname(target),tolerance=3e-12)
    expect_identical(dim(fit$vcov),c(ncol(xx),ncol(xx)))
    expect_identical(fit$vcov,t(fit$vcov))
    expect_true(min(eigen(fit$vcov,symmetric=TRUE,only.values=TRUE)$values)>-3e-12)
    expect_identical(fit$se,sqrt(diag(fit$vcov)))
  }
  fit <- solve.linear(y,X,y+yhat,yhat,1L)
  scaled <- solve.linear(y,sweep(X,2L,c(2,.5),`*`),y+yhat,yhat,1L)
  expect_equal(scaled$vcov,fit$vcov/outer(c(2,.5),c(2,.5)),tolerance=3e-12)
  expect_equal(scaled$coef,fit$coef/c(2,.5),tolerance=3e-12)
  expect_error(solve.linear(y,cbind(X[,1L],X[,1L]),y+yhat,yhat,1L),
               "residualized linear regressors are rank deficient",fixed=TRUE)
})

test_that("public partial-linear covariance matches the known-error witness", {
  if(exists("spawn_mpi_slaves",mode="function") && !spawn_mpi_slaves())
    skip("Could not initialize MPI context")
  g <- expand.grid(z=letters[1:2],r=c(-2,-1,1,2),e=c(-1,1),rep=1:2)
  X <- data.frame(x=g$r); Z <- data.frame(z=factor(g$z)); n <- nrow(X)
  for(hetero in c(FALSE,TRUE)) {
    u <- if(hetero) g$r*g$e else sqrt(2.5)*g$e
    y <- 2*g$r+ifelse(g$z=="a",1,3)+u
    bw <- npplregbw(xdat=X,ydat=y,zdat=Z,bws=matrix(0,2,1),bandwidth.compute=FALSE)
    fit <- npplreg(bws=bw,txdat=X,tydat=y,tzdat=Z,residuals=TRUE,se=TRUE)
    target <- n/(n-2)*sum(g$r^2*u^2)/sum(g$r^2)^2
    expect_equal(as.double(coef(fit)),2,tolerance=3e-12)
    expect_equal(as.double(residuals(fit)),u,tolerance=3e-12)
    expect_equal(as.double(vcov(fit)),target,tolerance=3e-12)
    expect_identical(coef(fit,se=TRUE),sqrt(diag(vcov(fit))))
    if(hetero) expect_equal(as.double(vcov(fit)),.045333333333333333,tolerance=3e-12)
    else expect_equal(as.double(vcov(fit)),sum(u^2)/(n-2)/sum(g$r^2),tolerance=3e-12)
  }
})
