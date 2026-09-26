# Literal finite-sample GNN regression oracle, independent of package kernels.
gnn_tree_literal <- function(x, y, k, degree, basis, method, kernel, mixed = NULL) {
  x <- as.matrix(x); n <- nrow(x); p <- ncol(x)
  terms <- as.matrix(expand.grid(rep(list(0:degree), p)))
  if (basis == "glp") terms <- terms[rowSums(terms) <= degree,,drop=FALSE]
  if (basis == "additive") terms <- terms[rowSums(terms > 0) <= 1,,drop=FALSE]
  B <- vapply(seq_len(nrow(terms)), function(a)
    apply(sweep(x, 2, terms[a,], "^"), 1, prod), numeric(n))
  K <- if (kernel == "uniform") function(u) .5*(abs(u)<1) else
    function(u) 3/(4*sqrt(5))*pmax(0,1-u*u/5)
  fit <- lev <- numeric(n)
  for (j in seq_len(n)) {
    donors <- if (method == "cv.ls") setdiff(seq_len(n),j) else seq_len(n)
    w <- rep(1,length(donors))
    for (l in seq_len(p)) {
      dist <- sort(abs(x[-j,l]-x[j,l]))
      h <- if (k[l] <= length(dist)) dist[k[l]] else max(dist)*k[l]/length(dist)
      w <- w*K((x[j,l]-x[donors,l])/h)
    }
    if (!is.null(mixed)) {
      w <- w*ifelse(mixed$u[donors] == mixed$u[j], .8, .1)
      w <- w*.3^abs(mixed$o[donors]-mixed$o[j])
    }
    D <- B[donors,,drop=FALSE]; G <- crossprod(D,w*D)
    fit[j] <- sum(B[j,]*solve(G,crossprod(D,w*y[donors])))
    if (method == "cv.aic")
      lev[j] <- w[match(j,donors)]*sum(B[j,]*solve(G,B[j,]))
  }
  mse <- mean((y-fit)^2)
  if (method == "cv.ls") mse else log(mse)+(1+mean(lev))/(1-mean(lev)-2/n)
}

test_that("explicit GNN LP trees reach only the qualified search context", {
  old <- options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  ns <- asNamespace("np")
  b <- npregbw(xdat=data.frame(x=seq(-1,1,length.out=20)),ydat=sin(1:20),
    bws=8,bwtype="generalized_nn",regtype="lp",degree=1L,
    ckertype="epanechnikov",bandwidth.compute=FALSE)
  choose <- get(".npregbw_tree_code",ns)
  for (mode in list(FALSE,TRUE,"auto")) for (context in c(FALSE,TRUE)) {
    options(np.tree=mode)
    expect_identical(choose(b,1L,0L,"lp",cv.context=context),
      get(if(isTRUE(mode)&&context)"DO_TREE_YES"else"DO_TREE_NO",ns))
  }
  options(np.tree=TRUE)
  for (type in c("adaptive_nn","generalized_nn")) for (method in c("cv.ks","cv.check")) {
    z<-b;z$type<-type;z$method<-method
    expect_identical(choose(z,1L,0L,"lp",cv.context=TRUE),get("DO_TREE_NO",ns))
  }
  z<-b;z$type<-"adaptive_nn"
  expect_identical(choose(z,1L,0L,"lp",cv.context=TRUE),get("DO_TREE_YES",ns))
})

test_that("GNN LP tree objectives agree with literal radii and leverage", {
  old <- options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  set.seed(971);n<-96L
  data <- data.frame(x1=runif(n,-1,1),x2=runif(n,-1,1),x3=runif(n,-1,1))
  data[2,]<-data[1,] # genuine ties, with a positive NN radius
  y <- sin(2*data$x1)+.2*data$x2+rnorm(n,sd=.1)
  cases <- data.frame(p=c(1,1,1,1,2,2,3,3),d=c(0,1,2,3,2,2,1,2),
                      basis=c("glp","glp","glp","glp","additive","tensor","glp","glp"))
  ns<-asNamespace("np");evaluate<-get(".npregbw_eval_only",ns)
  for (row in seq_len(nrow(cases))) for (method in c("cv.ls","cv.aic"))
    for (kernel in c("epanechnikov","uniform")) {
      q<-cases[row,];x<-data[,seq_len(q$p),drop=FALSE];k<-rep(60L,q$p)
      expected<-gnn_tree_literal(x,y,k,q$d,q$basis,method,kernel)
      for (bernstein in c(FALSE,TRUE)) {
        options(np.tree=FALSE)
        b<-npregbw(xdat=x,ydat=y,bws=k,bwtype="generalized_nn",bwmethod=method,
          regtype="lp",degree=rep(q$d,q$p),basis=q$basis,
          bernstein.basis=bernstein,ckertype=kernel,bandwidth.compute=FALSE)
        off<-evaluate(x,y,b)$objective
        options(np.tree=TRUE);on<-evaluate(x,y,b)$objective
        expect_true(abs(on-expected)<2e-9,info=paste(row,method,kernel,bernstein,on,expected))
        expect_true(abs(on-off)<2e-9,info=paste("tree parity",row,method,kernel,bernstein))
      }
    }
})
