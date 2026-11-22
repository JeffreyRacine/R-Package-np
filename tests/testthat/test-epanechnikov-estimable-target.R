test_that("identified endpoint fits and hats agree with reduced weighted least squares", {
  old <- options(np.messages=FALSE, np.tree=FALSE)
  on.exit(options(old), add=TRUE)
  i <- seq_len(136L)
  x <- data.frame(u=rep(seq(0,4,by=.25),each=8L), v=.5*((i*.7548776662466927) %% 1))
  y <- cos(x$u)+2*x$v+.3*((i*.5698402909980532) %% 1 >= .5)+.05*sin(13*i)
  grid <- data.frame(u=c(0,4),v=.25)
  for (cols in list(1:2,2:1)) {
    bw <- npregbw(xdat=x[cols],ydat=y,bws=c(12,40)[cols],degree=c(3L,1L)[cols],
                  regtype="lp",basis="glp",bernstein.basis=TRUE,bwtype="generalized_nn",
                  ckertype="epanechnikov",bandwidth.compute=FALSE)
    fit <- npreg(bw,txdat=x[cols],tydat=y,exdat=grid[cols],gradients=TRUE,se=FALSE)
    H <- npreghat(bw,txdat=x[cols],exdat=grid[cols],s=c(0L,0L))
    D <- npreghat(bw,txdat=x[cols],exdat=grid[cols],s=c(0L,1L)[cols])
    for(j in 1:2) {
      u <- x$u-grid$u[j]
      v <- x$v-grid$v[j]
      weight <- pmax(0,1-u^2/(5*sort(abs(u))[12L]^2)) *
                pmax(0,1-v^2/(5*sort(abs(v))[40L]^2))
      # Three u levels identify the pure quadratic and u-by-v quadratic
      # columns. Omitting the aliased u^3 leaves both target functionals intact.
      A <- cbind(1,u,u^2,v,u*v,u^2*v)
      ref <- qr.coef(qr(A*sqrt(weight)),diag(sqrt(weight)))
      expect_lt(max(abs(H[j,]-ref[1,]),abs(D[j,]-ref[4,])),1e-8)
      expect_lt(max(abs(fit$mean[j]-sum(ref[1,]*y)),
                    abs(fit$grad[j,which(cols==2L)]-sum(ref[4,]*y))),1e-8)
    }
  }
})
