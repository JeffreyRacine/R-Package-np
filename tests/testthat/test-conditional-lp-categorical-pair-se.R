test_that("LP categorical covariance differentiates the paired accepted maps", {
  set.seed(721)
  n <- 40L
  x <- seq(-1, 1, length.out = n)
  y <- sin(x) + rnorm(n, sd = .4)
  z <- cbind(1, x, x^2)
  d <- rbind(c(1, .1, .01), c(1, .1, .01))
  w <- cbind(dnorm(x)*rep(c(.8,.2), n/2), dnorm(x)*rep(c(.2,.8), n/2))
  r <- matrix(dnorm((y-.3)/.6)/.6, n, 2L)
  native <- function(w, r) .Call("C_np_conditional_lp_pair_se", z, d, w, r, NULL)
  got <- native(w, r)
  map <- function(mass, j) {
    a <- crossprod(z, z*(mass*w[,j]))
    drop(d[j,] %*% solve(a, crossprod(z, mass*w[,j]*r[,j])))
  }
  eps <- 1e-5
  influence <- vapply(seq_len(n), function(i) {
    up <- down <- rep(1, n)
    up[i] <- 1+eps; down[i] <- 1-eps
    (map(up,1)-map(up,2)-map(down,1)+map(down,2))/(2*eps)
  }, 0.)
  expect_equal(got[1, 2:3], vapply(1:2, function(j) map(rep(1,n),j),0.), tolerance=1e-12)
  expect_equal(got[1, 1], sqrt(n/(n-1)*sum((influence-mean(influence))^2)), tolerance=1e-8)
  expect_identical(native(w[,c(1,1)],r)[1,1], 0)
  expect_lt(native(w,matrix(2,n,2))[1,1], 1e-13)
  expect_error(.Call("C_np_conditional_lp_pair_se", z, d, w[-1,], r, NULL), "dimensions")
})

test_that("unrequested conditional contrast uncertainty never enters its producer", {
  x <- data.frame(x=seq(.05,.95,length.out=24L),g=factor(rep(0:1,12L)))
  y <- data.frame(y=sin(x$x)+.1*as.integer(x$g))
  b <- npcdensbw(xdat=x,ydat=y,bws=c(.3,.3,.25),
    bandwidth.compute=FALSE,regtype="ll")
  testthat::local_mocked_bindings(.np_conditional_lp_pair_se=function(...)
    stop("unrequested paired producer"), .package="np")
  expect_no_error(npcdens(bws=b,txdat=x,tydat=y,se=TRUE))
  expect_no_error(npcdens(bws=b,txdat=x,tydat=y,gradients=TRUE,se=FALSE))
  expect_no_error(npcdens(bws=b,txdat=x,tydat=y,gradients=TRUE,se=TRUE,
                         .np_conditional_cat_se_demand=FALSE))
})

test_that("conditional LP categorical errors honor demand without changing points", {
  set.seed(813)
  n <- 40L
  x <- data.frame(x=runif(n,-1,1),g=factor(rep(0:1,n/2)),
                  o=ordered(rep(1:4,n/4)))
  y <- data.frame(y=.4*x$x+.2*as.integer(x$g)+rnorm(n,sd=.5))
  ex <- x[c(8,24,33),]; ey <- data.frame(y=c(.2,.5,.8))
  for (type in c("fixed","generalized_nn","adaptive_nn")) for (cdf in c(FALSE,TRUE)) {
    b <- (if(cdf) npcdistbw else npcdensbw)(xdat=x,ydat=y,
      bws=c(if(type=="fixed") .55 else 30,if(type=="fixed") .6 else 32,.25,.3),
      bandwidth.compute=FALSE,bwtype=type,regtype="lp",degree=2L)
    f <- if(cdf) npcdist else npcdens
    a <- f(bws=b,txdat=x,tydat=y,exdat=ex,eydat=ey,gradients=TRUE,se=TRUE)
    no <- f(bws=b,txdat=x,tydat=y,exdat=ex,eydat=ey,gradients=TRUE,se=TRUE,
            .np_conditional_cat_se_demand=FALSE)
    one <- f(bws=b,txdat=x,tydat=y,exdat=ex,eydat=ey,gradients=TRUE,se=TRUE,
             .np_conditional_cat_se_demand=c(TRUE,FALSE))
    expect_true(all(is.finite(a$congerr[,2:3])))
    expect_identical(a$congrad,no$congrad)
    expect_identical(fitted(a),fitted(no))
    expect_identical(a$congerr[,1],no$congerr[,1])
    expect_identical(a$congerr[,2],one$congerr[,2])
    expect_true(all(is.na(one$congerr[,3])))
    if (cdf) {
      q <- npqreg(bws=b,txdat=x,tydat=y,exdat=ex,tau=.5,gradients=TRUE,se=TRUE)
      expect_true(all(is.finite(q$quantgerr[,2:3])))
    }
  }
})
