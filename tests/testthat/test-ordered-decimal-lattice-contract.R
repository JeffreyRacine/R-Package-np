r21_ordered_oracle <- function(s, lambda, kernel, op = "normal") {
  distance <- abs(outer(s, s, "-"))
  raw <- lambda^distance
  normal <- switch(kernel,
    liracine = raw,
    nliracine = raw*(1-lambda)/(1+lambda),
    wangvanryzin = ifelse(distance == 0, 1-lambda, .5*(1-lambda)*raw),
    racineliyan = raw/rowSums(raw))
  if (op == "normal") return(normal)
  if (op == "integral") {
    if (kernel == "racineliyan") return(t(apply(normal, 1L, cumsum)))
    grid <- if (kernel == "liracine") seq(min(s), max(s)) else -120:120
    K <- lambda^abs(outer(s, grid, "-"))
    if (kernel == "nliracine") K <- K*(1-lambda)/(1+lambda)
    if (kernel == "wangvanryzin")
      K <- ifelse(outer(s,grid,"=="), 1-lambda, .5*(1-lambda)*K)
    return(vapply(s, function(q) rowSums(K[,grid <= q,drop=FALSE]), numeric(length(s))))
  }
  if (kernel == "racineliyan") return(tcrossprod(normal))
  grid <- -120:120
  K <- lambda^abs(outer(s,grid,"-"))
  if (kernel == "nliracine") K <- K*(1-lambda)/(1+lambda)
  if (kernel == "wangvanryzin")
    K <- ifelse(outer(s,grid,"=="), 1-lambda, .5*(1-lambda)*K)
  tcrossprod(K)
}

test_that("ordered kernels and caches preserve decimal-offset integer gaps", {
  old <- options(np.messages = FALSE, np.largeh = FALSE, np.largelambda = FALSE)
  on.exit(options(old), add = TRUE)
  s <- c(0,1,4,7)
  ids <- rep(seq_along(s), c(2,4,3,3))
  for (offset in c(0, .1, .212, -8.788))
    for (kernel in c("liracine","nliracine","wangvanryzin","racineliyan")) {
      labels <- s+offset
      d <- data.frame(o = ordered(labels[ids], levels=labels))
      e <- data.frame(o = ordered(labels, levels=labels))
      for (op in c("normal","integral","convolution")) {
        if (kernel == "liracine" && op == "convolution") next
        expected <- r21_ordered_oracle(s,.35,kernel,op)[ids,,drop=FALSE]
        got <- npksum(txdat=d,exdat=e,bws=.35,okertype=kernel,operator=op,
                     return.kernel.weights=TRUE)
        expect_equal(unname(got$kw),expected,tolerance=2e-12)
        expect_equal(as.numeric(got$ksum),colSums(expected),tolerance=2e-12)
      }
      eps <- 1e-6
      expected <- (r21_ordered_oracle(s,.35+eps,kernel)-
                   r21_ordered_oracle(s,.35-eps,kernel))/(2*eps)
      got <- npksum(txdat=d,exdat=e,bws=.35,okertype=kernel,
        compute.score=TRUE,return.kernel.weights=TRUE,
        return.derivative.kernel.weights=TRUE)
      expect_equal(as.numeric(got$p.kw),as.numeric(expected[ids,,drop=FALSE]),
                   tolerance=2e-9)
    }
})

test_that("ordered validation admits integer lattices but not fractional distances", {
  dlev <- getFromNamespace("dlev", "npRmpi")
  for (s in list(.1+0:9, .212+c(0,1,4,7), -8.788+c(0,1,4,7))) {
    x <- ordered(s,levels=s)
    expect_equal(dlev(x), s, tolerance=1e-15)
  }
  for (s in list(c(0,.5,1),c(.1,1.10000001,4.1),
                c(0,2,1),c(0,Inf),c(0,.Machine$integer.max))) {
    x <- ordered(s,levels=s)
    expect_error(dlev(x),"integer distances")
  }
  expect_identical(dlev(ordered(c("low","high"),levels=c("low","high"))),c(1,2))
  bridge <- function(s) .Call("C_np_ordered_rly_matrix", s,s,.35,s,PACKAGE="npRmpi")
  for (offset in c(.1,.212,-8.788)) {
    s <- c(0,1,4,7)
    expect_equal(bridge(s+offset),t(r21_ordered_oracle(s,.35,"racineliyan")),
                 tolerance=2e-13)
  }
  expect_error(bridge(c(.1,1.10000001,4.1)),"integer distances")
})

test_that("ordered lattice conversion reaches bandwidth families without recoding", {
  old <- options(np.messages = FALSE)
  on.exit(options(old),add=TRUE)
  s <- .1+0:3
  d <- data.frame(o=ordered(rep(s,each=6),levels=s),
                  x=seq(.1,.9,length.out=24))
  y <- sin(d$x)+as.integer(d$o)
  constructors <- list(
    function() npregbw(xdat=d,ydat=y,bws=c(.2,.3),bandwidth.compute=FALSE),
    function() npudensbw(dat=d,bws=c(.2,.3),bandwidth.compute=FALSE),
    function() npudistbw(dat=d,bws=c(.2,.3),bandwidth.compute=FALSE),
    function() npcdensbw(xdat=d,ydat=y,bws=c(.2,.3,.4),bandwidth.compute=FALSE),
    function() npcdistbw(xdat=d,ydat=y,bws=c(.2,.3,.4),bandwidth.compute=FALSE))
  for (construct in constructors) expect_no_error(construct())
})

test_that("ordered translation preserves public fits and uncertainty", {
  old <- options(np.messages = FALSE, np.largeh = FALSE, np.largelambda = FALSE)
  on.exit(options(old), add = TRUE)
  for (family in c("regression", "density", "distribution", "condensity", "condistribution")) {
    fits <- lapply(c(0, .212), function(offset) {
      s <- c(0, 1, 4, 7) + offset
      d <- data.frame(o = ordered(rep(s, each = 6), levels = s),
                      x = seq(.1, .9, length.out = 24))
      y <- sin(d$x) + as.integer(d$o)
      if (family == "regression") {
        b <- npregbw(xdat = d, ydat = y, bws = c(.2, .3),
                     regtype = "ll", bandwidth.compute = FALSE)
        npreg(b, txdat = d, tydat = y, gradients = TRUE, se = TRUE)
      } else if (family == "density") {
        b <- npudensbw(dat = d, bws = c(.2, .3), bandwidth.compute = FALSE)
        npudens(b, tdat = d, se = TRUE)
      } else if (family == "distribution") {
        b <- npudistbw(dat = d, bws = c(.2, .3), bandwidth.compute = FALSE)
        npudist(b, tdat = d, se = TRUE)
      } else if (family == "condensity") {
        b <- npcdensbw(xdat = d, ydat = y, bws = c(.2, .3, .4),
                       bandwidth.compute = FALSE)
        npcdens(b, txdat = d, tydat = y, gradients = TRUE, se = TRUE)
      } else {
        b <- npcdistbw(xdat = d, ydat = y, bws = c(.2, .3, .4),
                       bandwidth.compute = FALSE)
        npcdist(b, txdat = d, tydat = y, gradients = TRUE, se = TRUE)
      }
    })
    expect_equal(fitted(fits[[1L]]), fitted(fits[[2L]]), tolerance = 1e-12)
    expect_equal(se(fits[[1L]]), se(fits[[2L]]), tolerance = 1e-12)
    if (family %in% c("regression", "condensity", "condistribution"))
      expect_equal(gradients(fits[[1L]]), gradients(fits[[2L]]), tolerance = 1e-12)
  }
})
