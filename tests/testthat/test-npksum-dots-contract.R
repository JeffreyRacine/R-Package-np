test_that("npksum refuses ignored extras without evaluating them", {
  x <- seq(.1,.9,length.out=9)
  for(nm in c("sub", "subs", "na.act", "leave.one.outt", "bandwith.divide", "nosuch")) {
    hit <- 0L
    bad <- as.call(c(list(quote(npksum)),list(txdat=quote(x),bws=.2),
      setNames(list(quote({hit <- hit+1L; TRUE})),nm)))
    expect_error(eval(bad), "unused npksum argument")
    expect_identical(hit,0L)
  }
  d <- data.frame(x=x,y=x^2); hit <- 0L
  expect_error(npksum(y~x,data=d,bws=.2,nosuch={hit<-hit+1L;1}), "unused npksum")
  expect_identical(hit,0L)
  leaf <- getFromNamespace("npksum.default","npRmpi")
  expect_error(leaf(bws=.2,txdat=x,nosuch=stop("forced")), "unused npksum")
})

test_that("npksum matching, formula ownership and internal controls survive", {
  d <- data.frame(x=seq(.1,.9,length.out=12)); d$y <- d$x^2
  args <- list(txdat=d$x,bws=.2,ckertype="gaussian",leave.one.out=TRUE)
  full <- do.call(npksum,args)
  names(args)[3:4] <- c("ckerty","leave.one.o")
  expect_identical(do.call(npksum,args)$ksum,full$ksum)
  counts <- new.env(parent=emptyenv()); counts$subset <- counts$na <- 0L
  omit <- function(x,...) {counts$na <- counts$na+1L; stats::na.omit(x)}
  form <- npksum(y~x,data=d,bws=.2,sub={counts$subset<-counts$subset+1L;x>.3},na.act=omit)
  explicit <- npksum(txdat=d$x[d$x>.3],tydat=d$y[d$x>.3],bws=.2)
  expect_identical(form$ksum,explicit$ksum)
  expect_identical(counts$subset,1L);expect_identical(counts$na,1L)
  one <- npksum(txdat=d$x,bws=1)
  expect_identical(npksum(txdat=d$x,bws=1L)$ksum,one$ksum)
  stored <- npudensbw(~x,data=d,bws=1,bwscaling=FALSE,bandwidth.compute=FALSE)
  expect_identical(npksum(txdat=d$x,bws=stored)$ksum,one$ksum)
  power <- getFromNamespace(".npksum_power12","npRmpi")(bws=stored,txdat=d$x)
  expect_equal(as.vector(power$ksum),as.vector(one$ksum),tolerance=0)
  oracle <- rowSums(outer(d$x,d$x,function(a,b) dnorm(a-b)))
  expect_equal(as.vector(one$ksum),oracle,tolerance=1e-14)
  expect_error(npksum(txdat=d$x,bws=.2,cker="gaussian"), "matches multiple")
})
