# The reference deletes the observation before deriving either NN radius.
# It is deliberately independent of the package's NN and objective helpers.
cdf_deleted_radius <- function(z, at, k) sort(abs(z-at), partial=k)[k]
cdf_deleted_ann <- function(z, k) vapply(seq_along(z), function(j)
  cdf_deleted_radius(z[-j],z[j],k), 0.0)
cdf_deleted_value <- function(X,Y,b,grid,full=FALSE) {
  pkg <- getNamespaceName(environment(npcdistbw))
  get('.npcdistbw_eval_only',asNamespace(pkg))(
    xdat=X,ydat=Y,bws=b,gydat=grid,do.full.integral=full)$objective
}

test_that('beta conditional CDF uses deleted NN geometry on both sides', {
  old <- options(np.messages=FALSE,np.tree=FALSE,np.extendednn=FALSE)
  on.exit(options(old),add=TRUE)
  x <- c(.08,.12,.21,.32,.43,.54,.65,.73,.86,.94)
  y <- c(.14,.31,.09,.61,.35,.83,.48,.74,.92,.56)
  X <- data.frame(x=x);Y <- data.frame(y=y);k <- 4L
  for(type in c('generalized_nn','adaptive_nn')) {
    b <- npcdistbw(xdat=X,ydat=Y,bws=c(k,k),bwtype=type,
      bandwidth.compute=FALSE,cxkertype='beta',cykertype='beta',
      cxkerbound='fixed',cxkerlb=0,cxkerub=1,
      cykerbound='fixed',cykerlb=0,cykerub=1)
    for(full in c(FALSE,TRUE)) {
      grid <- if(full)NULL else data.frame(y=c(.2,y[3],.8))
      oracle <- mean(vapply(seq_along(x),function(i) {
        xx <- x[-i]; yy <- y[-i]; gg <- if(full)y[-i] else grid$y
        hx <- if(type=='generalized_nn')cdf_deleted_radius(xx,x[i],k) else cdf_deleted_ann(xx,k)
        wx <- dbeta(xx,1+x[i]/hx^2,1+(1-x[i])/hx^2)
        fit <- vapply(gg,function(t) {
          hy <- if(type=='generalized_nn')cdf_deleted_radius(yy,t,k) else cdf_deleted_ann(yy,k)
          sum(wx*pbeta(t,1+yy/hy^2,1+(1-yy)/hy^2))/sum(wx)
        },0.0)
        mean(((y[i]<=gg)-fit)^2)
      },0.0))
      expect_equal(cdf_deleted_value(X,Y,b,grid,full),oracle,tolerance=1e-10)
    }
  }
})

test_that('deleted CDF tiles cover wide and categorical responses with signed LP rows', {
  old <- options(np.messages=FALSE,np.tree=FALSE,np.extendednn=FALSE)
  on.exit(options(old),add=TRUE)
  set.seed(8301)
  n <- 16L;X <- data.frame(x=runif(n,.08,.92))
  continuous <- data.frame(y1=runif(n,.08,.92),y2=runif(n,.08,.92),y3=runif(n,.08,.92))
  for(type in c('generalized_nn','adaptive_nn')) for(d in c(0L,3L)) {
    Y <- continuous[seq_len(d)]
    Y$category <- ordered(rep(1:3,length.out=n))
    grid <- Y[c(2,7,11),,drop=FALSE]
    args <- list(xdat=X,ydat=Y,bws=c(rep(7,d),.2,7),bwtype=type,
      regtype='lp',degree=2L,bandwidth.compute=FALSE,
      cxkertype='beta',cxkerbound='fixed',cxkerlb=0,cxkerub=1)
    if(d)args <- c(args,list(cykertype='beta',cykerorder=4L,
      cykerbound='fixed',cykerlb=rep(0,d),cykerub=rep(1,d)))
    b <- do.call(npcdistbw,args)
    oracle <- mean(vapply(seq_len(n),function(i) {
      pred <- fitted(npcdist(bws=b,txdat=X[-i,,drop=FALSE],tydat=Y[-i,,drop=FALSE],
        exdat=X[rep(i,nrow(grid)),,drop=FALSE],eydat=grid))
      indicator <- Reduce(`&`,lapply(seq_along(Y),function(j)
        as.numeric(Y[[j]][i])<=as.numeric(grid[[j]])))
      mean((indicator-pred)^2)
    },0.0))
    expect_equal(cdf_deleted_value(X,Y,b,grid),oracle,tolerance=1e-10)
  }
})
