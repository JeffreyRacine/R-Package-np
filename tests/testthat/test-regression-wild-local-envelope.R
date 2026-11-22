test_that("Mammen local envelopes preserve the observed SE and use generation SD", {
  old.options <- options(np.messages=FALSE,np.plot.progress=FALSE)
  on.exit(options(old.options),add=TRUE)
  if (!spawn_mpi_slaves(1L)) skip("MPI pool unavailable in this context")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  set.seed(103141)
  x <- data.frame(x=seq(-1,1,length.out=80))
  y <- sin(3*x$x)+(1+x$x^2)*rnorm(80,sd=.3)
  bw <- npregbw(xdat=x,ydat=y,bws=.3,regtype="ll",bandwidth.compute=FALSE)
  fit <- npreg(bws=bw,txdat=x,tydat=y)
  grid <- data.frame(x=seq(-.8,.8,length.out=7))
  state <- .np_bias_regression_prepare(bw,x,y,np_bias_control(),adapt=TRUE,law="mammen")
  eval <- .npreg_publish_plot_rows(npreg(bws=bw,txdat=x,tydat=y,exdat=grid,
    se=TRUE,.np.fit.inference=TRUE))
  seed <- .Random.seed
  run <- function(cal) .np_with_seed(103142L,
    .np_bias_inference(state,grid,"mean",alpha=c(.1,.05),B=199L,law="mammen",
      fit=eval,calibration=cal,resolution.warn=FALSE))
  ref <- run("absolute")
  local <- run("basic-local-envelope")
  expect_identical(local$se,ref$se)
  expect_identical(local$center,ref$center)
  expect_identical(.np_bias_generation_scale(state),.np_bias_student_prepare(state)$root.variance)
  expect_true(any(abs(local$generation.se-local$se)>1e-8))
  for(b in names(local$bands)) {
    limits <- local$calibration$limits[[b]]
    expect_identical(unname(limits[,1,]),-unname(limits[,2,]))
    expect_true(all(local$bands[[b]][,1,2]<=local$bands[[b]][,1,1]))
    expect_true(all(local$bands[[b]][,2,2]>=local$bands[[b]][,2,1]))
  }
  expect_match(local$root.units,"generation roots",fixed=TRUE)
  expect_match(local$calibration$local.reference[[1]]$units,"generation-SD",fixed=TRUE)
  expect_identical(.Random.seed,seed)
  expect_error(run("basic-envelope"),"requires inid")
  expect_identical(.Random.seed,seed)
  public <- suppressWarnings(plot(fit,center="bias-corrected",errors="bootstrap",
    bootstrap="wild",band="all",B=199L,neval=7L,output="data",
    random.seed=103142L,interval="basic-local-envelope", boot.control=np_boot_control(wild="mammen")))
  d <- as.data.frame(public,band="simultaneous")
  expect_true(all(is.finite(d$lower)) && all(d$lower<=d$upper))
  expect_identical(.Random.seed,seed)
})

test_that("wild local envelope symmetry covers every target and alpha", {
  x <- array(seq(-12,11),c(4,2,3))
  x[,1,] <- -abs(x[,1,]);x[,2,] <- abs(x[,2,])
  z <- .np_band_tail_symmetrize(list(pointwise=x,simultaneous=x,bonferroni=x))
  for(b in names(z)) {
    expect_identical(unname(z[[b]][,1,]),-unname(z[[b]][,2,]))
    expect_true(all(z[[b]][,1,]<=x[,1,]))
    expect_true(all(z[[b]][,2,]>=x[,2,]))
  }
  expect_error(.np_band_tail_symmetrize(list(x=matrix(0,2,2))),"invalid root limits")
})

test_that("wild local endpoints match an independent complete-array oracle", {
  if (!spawn_mpi_slaves(1L)) skip("MPI pool unavailable in this context")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
# Independent full-array empirical-root oracle; no production tail/local helper used.
oracle <- function(root,center,scale,map,alpha,symmetric=FALSE){
 stopifnot(is.matrix(root),all(is.finite(root)),all(is.finite(scale)),all(scale>=0),nrow(root)==nrow(map));B<-ncol(root);M<-nrow(root)
 out<-setNames(lapply(1:3,function(i)array(NA_real_,c(M,2L,length(alpha)))),c('pointwise','simultaneous','bonferroni'))
 for(target in unique(map$target)){ii<-which(map$target==target);rr<-root[ii,,drop=FALSE];ss<-scale[ii];m<-length(ii)
  calc<-function(R){sr<-t(apply(R,1L,sort));ar<-t(apply(abs(R),1L,sort));lo<-sort(apply(R,2L,min));hi<-sort(apply(R,2L,max));ab<-sort(apply(abs(R),2L,max));ans<-setNames(lapply(1:3,function(i)array(0,c(m,2L,length(alpha)))),names(out))
   for(a in seq_along(alpha)){aa<-alpha[a];k<-1+floor(B*aa/2);kb<-1+floor(B*aa/(2*m));ka<-ceiling(B*(1-aa));kab<-ceiling(B*(1-aa/m));low<-lo[k];high<-hi[B+1-k];q<-ab[ka]
    ans$pointwise[,1,a]<-pmin(sr[,k],-ar[,ka]);ans$pointwise[,2,a]<-pmax(sr[,B+1-k],ar[,ka])
    ans$simultaneous[,1,a]<-min(low,-q);ans$simultaneous[,2,a]<-max(high,q)
    ans$bonferroni[,1,a]<-pmin(sr[,kb],low,-ar[,kab],-q);ans$bonferroni[,2,a]<-pmax(sr[,B+1-kb],high,ar[,kab],q)
   };ans}
  raw<-calc(rr);loc<-calc(rr/ifelse(ss==0,1,ss))
  for(b in names(out))for(a in seq_along(alpha)){
   low<-pmin(raw[[b]][,1,a],loc[[b]][,1,a]*ss);high<-pmax(raw[[b]][,2,a],loc[[b]][,2,a]*ss)
   if(symmetric){rad<-pmax(-low,high);low<- -rad;high<-rad}
   out[[b]][ii,1,a]<-center[ii]-high;out[[b]][ii,2,a]<-center[ii]-low
   if(any(ss==0)){jj<-ii[ss==0];stopifnot(all(root[jj,,drop=FALSE]==0));out[[b]][jj,1,a]<-out[[b]][jj,2,a]<-center[jj]}
  }
 };out
}
  # A fixed row provider isolates calibration from estimator arithmetic. Its
  # observed and generation scales differ, so substituting either is detectable.
  n <- 9L; B <- 199L; alpha <- c(.2,.05,.1)
  grid <- data.frame(x=seq_len(31L)); targets <- c("mean","d1")
  state <- list(x=data.frame(x=seq_len(n)),y=sin(seq_len(n)),
                variance=seq(.2,1.8,length.out=n),
                local.inference=list(generation=exp(seq(-2,2,length.out=n))),
                settings=list())
  W <- outer(seq_len(31L),seq_len(n),function(i,j)sin(i*j/9)/(1+i/31))
  W[1L,] <- 0
  zero <- FALSE
  provider <- function(ii,map) {
    w <- W[ii[map$row],,drop=FALSE]
    w[map$target=="d1",] <- -w[map$target=="d1",]
    if(zero) w[] <- 0
    factor <- sweep(w,2L,sqrt(state$variance),"*")
    list(W=w,center=drop(w%*%state$y),original_center=rep(0,nrow(w)),
         factor=factor,se=sqrt(rowSums(factor^2)),
         main.regularization=data.frame(row=map$row,ridge=0))
  }
  inference <- .np_bias_inference
  scope <- new.env(parent=environment(inference))
  scope$.np_bias_capture_prepare <- function(...)list(rows=provider)
  scope$.np_bias_capture_validate <- function(...)TRUE
  environment(inference) <- scope
  map <- expand.grid(row=seq_len(nrow(grid)),target=targets,stringsAsFactors=FALSE)
  value <- provider(seq_len(nrow(grid)),map)
  factor <- sweep(value$W,2L,sqrt(state$local.inference$generation),"*")
  generation.sd <- sqrt(rowSums(factor^2))
  # Independent full-array innovations and order statistics; no production
  # draw/tail/local reducer supplies the expected endpoints.
  innovations <- .np_with_seed(103162L, {
    u <- matrix(runif(n*B),n)
    ifelse(u<(sqrt(5)+1)/(2*sqrt(5)),(1-sqrt(5))/2,(1+sqrt(5))/2)
  })
  roots <- factor%*%innovations
  expected <- oracle(roots,value$center,generation.sd,map,alpha,TRUE)
  run <- function(law="mammen",s=state) .np_with_seed(103162L,
    inference(s,grid,targets,B=B,alpha=alpha,law=law,
              calibration="basic-local-envelope",resolution.warn=FALSE))
  set.seed(103161L);held <- .Random.seed
  result <- run()
  for(band in names(expected))
    expect_equal(unname(result$bands[[band]]),unname(expected[[band]]),tolerance=1e-12)
  expect_identical(result$se,value$se)
  expect_identical(result$generation.se,generation.sd)
  expect_equal(result$center,value$center,tolerance=1e-12)
  expect_identical(.Random.seed,held)
  # These controls establish that this fixture can detect a wrong scale or
  # omitted symmetry; it is not merely another finiteness/symmetry assertion.
  wrong.factor <- sweep(value$W,2L,sqrt(state$variance),"*")
  alternatives <- list(
    oracle(wrong.factor%*%innovations,value$center,sqrt(rowSums(wrong.factor^2)),map,alpha,TRUE),
    oracle(roots,value$center,value$se,map,alpha,TRUE),
    oracle(roots,value$center,generation.sd,map,alpha,FALSE))
  for(wrong in alternatives)
    expect_gt(max(abs(unlist(expected)-unlist(wrong))),1e-5)
  for(law in c("rademacher","gaussian")) {
    expect_error(run(law=law),"requires inid")
    expect_identical(.Random.seed,held)
  }
  categorical <- state;categorical$categorical.only <- TRUE
  expect_error(run(s=categorical),"categorical-only")
  expect_identical(.Random.seed,held)
  zero <- TRUE
  empty <- run()
  expect_true(all(unlist(empty$bands)==0))
  expect_true(all(empty$se==0))
  expect_true(all(empty$generation.se==0))
  expect_identical(.Random.seed,held)
})


test_that("public wild local admission keeps its law and categorical exclusions", {
  if (!spawn_mpi_slaves(1L)) skip("MPI pool unavailable in this context")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old.options <- options(np.messages=FALSE,np.plot.progress=FALSE)
  on.exit(options(old.options),add=TRUE)
  set.seed(103163L)
  x <- data.frame(x=seq(-1,1,length.out=40))
  y <- sin(3*x$x)+cos(seq_len(40))*.2
  bw <- npregbw(xdat=x,ydat=y,bws=.3,regtype="ll",bandwidth.compute=FALSE)
  fit <- npreg(bws=bw,txdat=x,tydat=y)
  categorical <- data.frame(u=factor(rep(c("a","b"),20)))
  cbw <- npregbw(xdat=categorical,ydat=y,bws=.2,regtype="lc",bandwidth.compute=FALSE)
  cfit <- npreg(bws=cbw,txdat=categorical,tydat=y)
  held <- .Random.seed
  run <- function(object,calibration,wild) suppressWarnings(plot(object,
    errors="bootstrap",center="bias-corrected",bootstrap="wild",B=99L,
    neval=5L,band="all",output="data",random.seed=103164L,
    interval=regression_test_interval(calibration), boot.control=np_boot_control(wild=wild)))
  expect_error(run(fit,"basic-local-envelope","rademacher"),"requires inid")
  expect_identical(.Random.seed,held)
  expect_error(run(cfit,"basic-local-envelope","mammen"),"categorical-only")
  expect_identical(.Random.seed,held)
  incumbent <- run(cfit,"symmetric-sd","mammen")
  d <- as.data.frame(incumbent,band="simultaneous")
  expect_true(all(is.finite(d$lower)) && all(d$lower<=d$upper))
  expect_identical(.Random.seed,held)
})
