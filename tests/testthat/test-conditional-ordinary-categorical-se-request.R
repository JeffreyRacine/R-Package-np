test_that("ordinary categorical errors use complete paired endpoint influences", {
  i <- seq_len(25L)
  x <- data.frame(o=ordered(rep(c("a","b","c"),length.out=25L)),
                  x=.04+.92*i/26,
                  u=factor(rep(c("b","c","a","a"),length.out=25L),levels=c("a","b","c")))
  y <- data.frame(y=.03+.94*((i*7L)%%29)/29)
  ex <- x[c(1L,8L,19L),,drop=FALSE]
  ex$o <- ordered(c("a","b","c"),levels=levels(x$o))
  ex$u <- factor(c("a","b","c"),levels=levels(x$u))
  ey <- data.frame(y=c(.25,.47,.73))
  b <- npcdensbw(xdat=x,ydat=y,bws=c(.16,.25,.3,.2),
    bandwidth.compute=FALSE,bwscaling=FALSE,regtype="lc",
    oxkertype="wangvanryzin",uxkertype="aitchisonaitken")
  run <- function(demand=NULL) npcdens(bws=b,txdat=x,tydat=y,
    exdat=ex,eydat=ey,gradients=TRUE,.np_conditional_cat_se_demand=demand)
  fit <- run(); none <- run(FALSE); subset <- run(c(TRUE,FALSE))
  expected <- matrix(0,3L,2L)
  for(j in seq_len(3L)) {
    cx <- dnorm((ex$x[j]-x$x)/.3)/.3
    response <- dnorm((ey$y[j]-y$y)/.16)/.16
    endpoint <- function(o,u) {
      delta <- abs(as.integer(x$o)-o)
      D <- cx*ifelse(delta==0,.75,.375*.25^delta)*
        ifelse(as.integer(x$u)==u,.8,.1)
      N <- D*response; m <- sum(N)/sum(D)
      (N-m*D)/sum(D)
    }
    o <- as.integer(ex$o[j]); u <- as.integer(ex$u[j])
    a <- endpoint(o,u)
    expected[j,1L] <- sqrt(25/24*sum((a-endpoint(if(o==1L)2L else o-1L,u))^2))
    expected[j,2L] <- sqrt(25/24*sum((a-endpoint(o,1L))^2))
  }
  expect_equal(unname(fit$congerr[,c(1L,3L)]),expected,tolerance=3e-10)
  expect_identical(fitted(fit),fitted(none))
  expect_identical(se(fit),se(none))
  expect_identical(gradients(fit),gradients(none))
  expect_identical(subset$congerr[,1L],fit$congerr[,1L])
  expect_identical(subset$congerr[,3L],none$congerr[,3L])
  expect_identical(fit$congerr[,2L],none$congerr[,2L])
  expect_identical(unname(fit$congerr[1L,3L]),0)
  expect_error(run(c(TRUE,NA)),"demand")
  expect_error(npcdens(bws=b,txdat=x,tydat=y,gradients=TRUE,
    .np_conditional_cat_se_demnad=FALSE),"unrecognized|unused|unknown")
})

test_that("ordinary categorical inference retains legitimate zero variance", {
  x <- data.frame(u=factor(rep(c("a","b"),10L)))
  y <- data.frame(y=seq(.1,.9,length.out=20L))
  b <- npcdensbw(xdat=x,ydat=y,bws=c(.2,.2),bandwidth.compute=FALSE,
    regtype="lc",cykertype="uniform")
  fit <- npcdens(bws=b,txdat=x,tydat=y,exdat=x[1:2,,drop=FALSE],
    eydat=data.frame(y=c(10,10)),gradients=TRUE)
  expect_identical(as.vector(fitted(fit)),c(0,0))
  expect_identical(as.vector(fit$congerr),c(0,0))
})
