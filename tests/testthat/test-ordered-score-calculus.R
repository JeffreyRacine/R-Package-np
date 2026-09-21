ordered_score_oracle <- function(distance, lambda, kernel) {
  dnum <- ifelse(distance == 0, 0,
                 distance * lambda^pmax(0, distance - 1))
  if (kernel == "wangvanryzin")
    return(ifelse(distance == 0, -1,
      0.5 * lambda^pmax(0, distance - 1) *
        (distance - (distance + 1) * lambda)))
  dnum * (1-lambda)/(1+lambda) -
    2 * lambda^distance/(1+lambda)^2
}

test_that("ordered bandwidth scores differentiate their normal kernels", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages=FALSE, np.largelambda=FALSE)
  on.exit(options(old), add=TRUE)
  for (values in list(0:3, c(0, 2, 5))) {
    train <- rep(values, each=2L)
    x <- data.frame(o=ordered(train, levels=values))
    e <- data.frame(o=ordered(values, levels=values))
    distance <- abs(outer(train, values, "-"))
    for (kernel in c("wangvanryzin", "nliracine")) {
      for (lambda in c(0, .01, .42, .99, 1)) {
        actual <- npksum(txdat=x, exdat=e, bws=lambda,
                          okertype=kernel, compute.score=TRUE)
        expected <- colSums(ordered_score_oracle(distance,lambda,kernel))
        expect_equal(as.double(actual$p.ksum), expected, tolerance=2e-13)
        ordinary <- npksum(txdat=x, exdat=e, bws=lambda, okertype=kernel)
        expect_equal(actual$ksum, ordinary$ksum, tolerance=0)
      }
      delta <- 1e-6
      plus <- npksum(txdat=x, exdat=e, bws=.42+delta, okertype=kernel)$ksum
      minus <- npksum(txdat=x, exdat=e, bws=.42-delta, okertype=kernel)$ksum
      expected <- colSums(ordered_score_oracle(distance,.42,kernel))
      expect_equal(as.double((plus-minus)/(2*delta)), expected, tolerance=1e-8)
    }
  }
})

test_that("ordered scores preserve selected-coordinate packing in mixed data", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages=FALSE,np.largelambda=FALSE,np.largeh=FALSE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(a=ordered(c(0,1,2,3,0,2),levels=0:3),
                  z=c(-.6,-.3,.1,.4,.8,1),
                  b=ordered(c(3,2,1,0,2,1),levels=0:3))
  e <- x[c(1,4,6),,drop=FALSE]
  normal <- function(d,h,k) if(k=="wangvanryzin")
    ifelse(d==0,1-h,.5*(1-h)*h^d) else h^d*(1-h)/(1+h)
  d1 <- abs(outer(as.numeric(as.character(x$a)),
                  as.numeric(as.character(e$a)),"-"))
  d2 <- abs(outer(as.numeric(as.character(x$b)),
                  as.numeric(as.character(e$b)),"-"))
  continuous <- dnorm(outer(x$z,e$z,"-")/.4)
  for(k in c("wangvanryzin","nliracine")) {
    actual <- npksum(txdat=x,exdat=e,bws=c(.37,.4,.61),
                      okertype=k,compute.score=TRUE)
    expected <- cbind(
      colSums(ordered_score_oracle(d1,.37,k)*continuous*normal(d2,.61,k)),
      colSums(normal(d1,.37,k)*continuous*ordered_score_oracle(d2,.61,k)))
    expect_equal(unname(as.matrix(actual$p.ksum)),expected,tolerance=2e-13)
  }
})
