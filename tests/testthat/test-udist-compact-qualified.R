test_that("Epanechnikov GNN CDF includes the variation in each query radius", {
  x <- seq(-1,1,length.out=53)^3 + seq(0,.18,length.out=53)
  dat <- data.frame(x);g <- data.frame(x=c(-.49,.04,.55))
  b <- npudistbw(dat=dat,bwtype="generalized_nn",ckertype="epanechnikov",
                 bws=12,bandwidth.compute=FALSE)
  st <- .np_udist_bc_prepare(b,dat,allow.compact=TRUE)
  rows <- .np_udist_bc_rows(st,g)
  dist <- abs(outer(g$x,x,"-"));k<-st$rank$lookup
  w<-min(floor(length(x)^(2/3)),floor(min(k,length(x)-k)/2))
  order.stats<-t(apply(dist,1,sort))
  r<-order.stats[,k];den<-(2*w/length(x))/(order.stats[,k+w]-order.stats[,k-w])
  H <- function(t) { t[] <- pmax(-1,pmin(1,t)); .5+(45*t-50*t^3+21*t^5)/32 }
  A <- outer(g$x,x,"-")/sqrt(5);J <- H(A/r)
  eps <- 1e-6;deriv <- (rowMeans(H(A/(r+eps)))-rowMeans(H(A/(r-eps))))/(2*eps)
  phi <- J-(deriv/den)*(dist<=r)
  F <- (phi-rowMeans(phi))/length(x)
  expect_equal(rows$center,rowMeans(J),tolerance=1e-12)
  expect_lt(max(abs(rows$factor-F)),2e-9)
  expect_gt(max(abs(F-(J-rowMeans(J))/length(x))),1e-4)
  expect_error(.np_udist_joint_inference(st,g,law="inid",B=19L),"asymptotic calibration")
  wrong <- b;wrong$ckertype<-"uniform"
  expect_error(.np_udist_bc_prepare(wrong,dat,allow.compact=TRUE),"Gaussian")
})

