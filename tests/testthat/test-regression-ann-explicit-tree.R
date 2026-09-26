# Independent donor-specific delete-one radii and weighted LP oracle.
ann_tree_literal <- function(x,y,k,d,method,kernel,basis="glp",mixed=NULL){
 x<-as.matrix(x);n<-nrow(x);p<-ncol(x);terms<-as.matrix(expand.grid(rep(list(0:d),p)))
 if(basis=="glp")terms<-terms[rowSums(terms)<=d,,drop=FALSE]
 if(basis=="additive")terms<-terms[rowSums(terms>0)<=1,,drop=FALSE]
 B<-vapply(seq_len(nrow(terms)),function(a)apply(sweep(x,2,terms[a,],"^"),1,prod),numeric(n))
 K<-if(kernel=="uniform")function(u).5*(abs(u)<1)else function(u)3/(4*sqrt(5))*pmax(0,1-u*u/5)
 fit<-lev<-numeric(n)
 for(j in seq_len(n)){
  donors<-if(method=="cv.ls")setdiff(seq_len(n),j)else seq_len(n);w<-rep(1,length(donors))
  for(q in seq_along(donors)){i<-donors[q];for(l in seq_len(p)){
   dist<-sort(abs(x[setdiff(donors,i),l]-x[i,l]));h<-if(k[l]<=length(dist))dist[k[l]]else max(dist)*k[l]/length(dist)
   w[q]<-w[q]*K((x[j,l]-x[i,l])/h)/h
  }}
  if(!is.null(mixed))w<-w*ifelse(mixed$u[donors]==mixed$u[j],.8,.1)*.3^abs(mixed$o[donors]-mixed$o[j])
  D<-B[donors,,drop=FALSE];G<-crossprod(D,w*D);fit[j]<-sum(B[j,]*solve(G,crossprod(D,w*y[donors])))
  if(method=="cv.aic")lev[j]<-w[match(j,donors)]*sum(B[j,]*solve(G,B[j,]))
 }
 mse<-mean((y-fit)^2);if(method=="cv.ls")mse else log(mse)+(1+mean(lev))/(1-mean(lev)-2/n)
}

test_that("ANN LP tree moments preserve ordinary CVLS and CVAIC objectives", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(552)
  n <- 80L
  dat <- data.frame(x = runif(n, -1, 1), z = runif(n, -1, 1))
  dat[2, ] <- dat[1, ] # positive radii despite a genuine tie
  y <- sin(dat$x) + .3 * dat$z + rnorm(n, sd = .1)
  evaluate <- getFromNamespace(".npregbw_eval_only", "np")
  for (p in 1:2) for (method in c("cv.ls", "cv.aic"))
    for (kernel in c("epanechnikov", "uniform")) {
      x <- dat[seq_len(p)]
      k <- rep(60L, p)
      expected <- ann_tree_literal(x, y, k, p, method, kernel)
      for (bernstein in c(FALSE, TRUE)) {
        options(np.tree = FALSE)
        b <- npregbw(xdat = x, ydat = y, bws = k,
          bwtype = "adaptive_nn", bwmethod = method, regtype = "lp",
          degree = rep(p, p), bernstein.basis = bernstein,
          ckertype = kernel, bandwidth.compute = FALSE)
        off <- evaluate(x, y, b)$objective
        options(np.tree = TRUE)
        on <- evaluate(x, y, b)$objective
        info <- paste(p, method, kernel, bernstein)
        expect_true(abs(on - expected) < 2e-9, info = info)
        expect_true(abs(on - off) < 2e-9, info = info)
      }
    }
})

