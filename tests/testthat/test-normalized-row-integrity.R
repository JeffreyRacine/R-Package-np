.n12_pkg <- function() if ("package:npRmpi" %in% search()) "npRmpi" else "np"
.n12_owner <- function(name) getFromNamespace(name, .n12_pkg())
.n12_local <- function(expr) {
  if (.n12_pkg() == "npRmpi")
    .n12_owner(".npRmpi_with_local_regression")(expr)
  else force(expr)
}
.n12_capture <- function(expr) {
  notices <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    notices <<- c(notices, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, notices = notices)
}

test_that("normalization retains finite signed and representable small sums", .n12_local({
  den <- c(2, 2^-70, -2^-70, 2^-1072)
  got <- .n12_owner(".np_normalization_denominator")(den, "test")
  expect_identical(got, den)
  H <- rbind(c(.25,.75), c(-.5,1.5), c(-.25,-.75)) * 2^-60
  y <- c(3,7)
  counts <- matrix(c(2,0, 1,1, 0,2), nrow=2)
  den <- crossprod(counts,t(H))
  expected <- crossprod(counts,t(H)*y) / den
  boot <- .n12_owner(".np_inid_lc_boot_from_hat")(
    H=H, ydat=y, B=3L, counts=counts)
  expect_identical(boot$t, expected)
  expect_identical(boot$t0, as.vector(H %*% y))
  ordinary <- .n12_owner(".np_inid_lc_boot_from_hat")(
    H=H / 2^-60, ydat=y, B=3L, counts=counts)
  expect_equal(boot$t, ordinary$t, tolerance=1e-14)
  expect_identical(dim(boot$t), c(3L,3L))
}))

test_that("required ratios fail without dropping or replacing a bootstrap draw", .n12_local({
  normalize <- .n12_owner(".np_normalization_denominator")
  for (d in list(c(1,0), c(1,Inf), c(1,NA_real_))) {
    expect_error(normalize(d,"required"), "normalizing weight sum")
  }
  expect_error(normalize(c(1,0),"signed",TRUE,zero.rows=c(FALSE,FALSE)),
               "zero normalizing weight sum")
  allowed <- normalize(c(1,0),"external",TRUE,zero.rows=c(FALSE,TRUE))
  expect_identical(allowed,c(1,NA_real_))
  set.seed(731)
  before <- .Random.seed
  H <- matrix(c(1,0), nrow=1)
  counts <- matrix(c(1,1, 0,2), nrow=2)
  expect_error(.n12_owner(".np_inid_lc_boot_from_hat")(
    H=H,ydat=c(4,7),B=2L,counts=counts),
    "bootstrap replication 2, evaluation row 1")
  expect_identical(.Random.seed,before)
  expect_error(.n12_owner(".np_bootstrap_ratio")(
    matrix(1,2,1),matrix(c(1,0),2,1),"chunk",first.replication=5L),
    "bootstrap replication 6, evaluation row 1")
}))

test_that("public index hats keep constant reproduction and external row identity", .n12_local({
  old <- options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(x=seq(0,1,length.out=32L), z=sin(seq_len(32L)))
  y <- sin(x$x*3)
  for (type in c("fixed","generalized_nn","adaptive_nn")) {
    h <- if(type=="fixed") .15 else 4
    for (spec in list(list(regtype="lc"), list(regtype="lp",degree=0L))) {
      b <- do.call(npindexbw,c(list(xdat=x,ydat=y,bws=c(1,0,h),
        bandwidth.compute=FALSE,bwtype=type,ckertype="gaussian"),spec))
      ex <- data.frame(x=c(.4,1.5,3),z=0)
      H <- npindexhat(b,txdat=x,exdat=ex)
      a <- npindexhat(b,txdat=x,exdat=ex,y=rep(7,32),output="apply")
      expect_equal(rowSums(H),rep(1,3),tolerance=1e-12)
      expect_equal(a,rep(7,3),tolerance=1e-12)
      expect_equal(as.vector(H %*% rep(7,32)),a,tolerance=1e-12)
      expect_null(attr(a,".np.empty.rows",exact=TRUE))
    }
  }
  b <- npindexbw(xdat=x,ydat=y,bws=c(1,0,.15),bandwidth.compute=FALSE,
                  regtype="lc",ckertype="gaussian")
  for (order in c(4L,6L,8L)) {
    signed <- npindexbw(xdat=x,ydat=y,bws=c(1,0,.3),bandwidth.compute=FALSE,
                        regtype="lc",ckertype="gaussian",ckerorder=order)
    sx <- data.frame(x=c(.4,3),z=0)
    expect_equal(rowSums(npindexhat(signed,txdat=x,exdat=sx)),c(1,1),
                 tolerance=1e-12)
    expect_equal(npindexhat(signed,txdat=x,exdat=sx,y=rep(7,32),output="apply"),
                 c(7,7),tolerance=1e-12)
  }
  ex <- data.frame(x=c(.4,100,.6),z=0)
  H <- .n12_capture(npindexhat(b,txdat=x,exdat=ex))
  a <- .n12_capture(npindexhat(b,txdat=x,exdat=ex,y=rep(7,32),output="apply"))
  expect_length(H$notices,1L)
  expect_length(a$notices,1L)
  expect_match(a$notices,"all computed kernel weights")
  expect_false(grepl("outside.*support|underflow|ties",a$notices))
  expect_true(all(is.na(H$value[2,])))
  expect_identical(is.na(a$value),c(FALSE,TRUE,FALSE))
  expect_equal(a$value[c(1,3)],c(7,7),tolerance=1e-12)
  expect_null(attr(H$value,".np.empty.rows",exact=TRUE))
  constraint <- .n12_capture(npindexhat(b,txdat=x,exdat=ex,y=rep(7,32),
                                       output="constraint"))
  expect_length(constraint$notices,1L)
  expect_equal(constraint$value,t(H$value)*7,tolerance=1e-12)
  index.train <- data.frame(index=x$x)
  index.eval <- data.frame(index=ex$x)
  expect_error(.n12_owner(".np_indexhat_exact")(
    b,index.train,index.eval,y=rep(7,32),output="apply"),
    "normalizing weight sum")
  point <- .n12_owner(".np_plot_singleindex_hat_apply_index")(
    b,index.train,index.eval,rep(7,32),allow.empty.rows=TRUE)
  expect_identical(attr(point,".np.empty.rows",exact=TRUE),c(0L,1L,0L))
}))

test_that("conditional ratio empty X and zero Y have different contracts", .n12_local({
  old <- options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(x=seq(0,1,length.out=32L))
  y <- data.frame(y=sin(x$x*3))
  b <- npcdensbw(xdat=x,ydat=y,bws=c(4,4),bwtype="adaptive_nn",
    bandwidth.compute=FALSE,regtype="lc",nomad=FALSE,cxkertype="uniform")
  ex <- data.frame(x=c(.4,2))
  H <- .n12_capture(npcdenshat(b,txdat=x,tydat=y,exdat=ex,eydat=c(.5,.5)))
  a <- .n12_capture(npcdenshat(b,txdat=x,tydat=y,exdat=ex,eydat=c(.5,.5),
                              y=rep(1,32),output="apply"))
  expect_equal(sum(grepl("all computed kernel weights",H$notices)),1L)
  expect_equal(sum(grepl("all computed kernel weights",a$notices)),1L)
  expect_true(all(is.na(H$value[2,])))
  expect_equal(a$value[1],sum(H$value[1,]),tolerance=1e-12)
  expect_true(is.na(a$value[2]))
  expect_null(attr(a$value,".np.empty.rows",exact=TRUE))
  zero.y <- npcdensbw(xdat=x,ydat=y,bws=c(4,4),bwtype="adaptive_nn",
    bandwidth.compute=FALSE,regtype="lc",nomad=FALSE,cykertype="uniform")
  z <- suppressWarnings(npcdenshat(zero.y,txdat=x,tydat=y,
    exdat=data.frame(x=.4),eydat=100,y=rep(1,32),output="apply"))
  expect_identical(z,0)
  proper <- .n12_owner(".npConmodeProperProbabilities")(
    rbind(c(.4,.6),c(NA_real_,NA_real_)),levels=c("a","b"))
  expect_true(all(is.na(proper$probabilities[2,])))
  expect_equal(proper$proper.info$invalid.rows,1L)
  expect_identical(proper$probabilities[1,],c(a=.4,b=.6))
}))

test_that("signed moment cancellation is not mislabeled an empty external row", .n12_local({
  pkg <- .n12_pkg()
  run <- function() .n12_owner(".np_indexhat_lc_moment_apply")(
    matrix(c(1,2),ncol=1), list(txdat=data.frame(x=c(0,1)),
      exdat=data.frame(x=2),ckertype="gaussian",ckerorder=4L,bwtype="fixed"),
    allow.empty.rows=TRUE)
  probe <- function(...) list(kw=matrix(c(1,-1),ncol=1))
  testthat::with_mocked_bindings({
    expect_error(run(),"zero normalizing weight sum")
  }, .np_index_kernel_moments=function(...) list(numerator=matrix(1,1,1),
                                                denominator=0),
     .np_index_kernel_sum=probe, .package=pkg)
}))
