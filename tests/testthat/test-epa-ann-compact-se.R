test_that("sparse Epanechnikov training SEs retain the fitted points and full hat variance", {
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force=TRUE), add=TRUE)
  v <- readRDS(test_path("fixtures","r112-epanechnikov-ann.rds"))
  x <- data.frame(x=v$x)
  bw <- npregbw(xdat=x, ydat=v$y, bws=53, bandwidth.compute=FALSE,
    regtype="lp", degree=3, basis="glp", bernstein.basis=TRUE,
    bwtype="adaptive_nn", ckertype="epanechnikov", ckerorder=2)
  point <- npreg(bw,txdat=x,tydat=v$y,gradients=TRUE)
  fit <- npreg(bw,txdat=x,tydat=v$y,gradients=TRUE,se=TRUE)
  expect_identical(fit$mean,point$mean)
  expect_identical(fit$grad,point$grad)
  H <- npreghat(bw,txdat=x,exdat=x,s=0L)
  D <- diag(nrow(x))-H
  z <- drop(D %*% v$y)/sqrt(rowSums(D^2))
  take <- c(which.min(v$x),which.max(v$x))
  G <- npreghat(bw,txdat=x,exdat=x[take,,drop=FALSE],s=1L)
  expect_lt(max(abs(fit$merr[take]-sqrt(drop(H[take,]^2 %*% (z^2))))),1e-8)
  expect_lt(max(abs(fit$gerr[take,1]-sqrt(drop(G^2 %*% (z^2))))),1e-8)
})
