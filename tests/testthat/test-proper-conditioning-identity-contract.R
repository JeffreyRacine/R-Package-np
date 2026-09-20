test_that("conditioning identity does not depend on coordinate origin or labels", {
  for (origin in c(0,1e9)) {
    x <- data.frame(x=rep(origin+c(.2,.25),each=3L))
    expect_identical(as.integer(.np_condens_make_x_groups(x)),rep(1:2,each=3L))
    for (group in list(.np_condens_slice_groups,.np_condist_slice_groups))
      expect_equal(unname(group(x)),list(1:3,4:6))
  }
  # Distinct representable doubles may have the same formatted factor label.
  x <- data.frame(x=rep(c(1e9,1e9+2^-22),each=2L),
                  f=factor(rep(c("b","a"),2L)))
  expect_identical(as.integer(.np_condens_make_x_groups(x)),1:4)
  for (group in list(.np_condens_slice_groups,.np_condist_slice_groups)) {
    out <- group(x)
    expect_equal(length(out),4L)
    expect_equal(sort(unlist(out,use.names=FALSE)),1:4)
  }
  # Both selectors retain their ordinary repeated-row behavior.
  x <- data.frame(x=c(.2,.1,.2,.1),f=factor(c("a","b","a","b")))
  expect_equal(unname(.np_condens_slice_groups(x)),list(c(2L,4L),c(1L,3L)))
})

test_that("translated conditional properization still has two fixed-X slices", {
  old <- options(np.messages=FALSE)
  on.exit(options(old),add=TRUE)
  set.seed(920326)
  x <- data.frame(x=runif(60,.1,.4))
  y <- data.frame(y=10*x$x+rnorm(60,sd=.3))
  e <- data.frame(x=rep(c(.2,.25),each=21L))
  ey <- data.frame(y=rep(seq(-1,6,length.out=21L),2L))
  for (family in c("npcdens","npcdist")) for (origin in c(0,1e9)) {
    ctrl <- list(fail.on.unsupported=TRUE)
    if (family=="npcdens") ctrl$mass.warn.tol <- 0
    fit <- do.call(get(family),list(txdat=x+origin,tydat=y,exdat=e+origin,
      eydat=ey,bws=c(.3,.055),regtype="lp",degree=1L,bernstein.basis=TRUE,
      se=FALSE,proper=TRUE,proper.control=ctrl))
    expect_true(fit$proper.applied)
    expect_equal(fit$proper.info$slice.count,2L)
    expect_true(all(is.finite(fitted(fit))))
  }
})
