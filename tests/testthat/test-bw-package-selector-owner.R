test_that("package selector heads ignore caller masking without moving arguments", {
  ns <- environment(npudens)
  package <- getNamespaceName(ns)
  selectors <- c("npregbw", "npplregbw", "npudensbw", "npudistbw", "npcdensbw",
                 "npcdistbw", "npindexbw", "npscoefbw")
  caller <- new.env(parent = ns)
  hits <- 0L
  for (name in selectors)
    makeActiveBinding(name, function(value) {
      hits <<- hits + 1L; stop("caller selector must not be read")
    }, caller)
  for (name in selectors) {
    call <- as.call(list(as.name(name), xdat = quote(stop("do not force data"))))
    owned <- .np_bw_package_selector_call(call)
    expect_identical(owned[[1L]], substitute(PKG::FN,
      list(PKG = as.name(package), FN = as.name(name))))
    expect_identical(as.list(owned)[-1L], as.list(call)[-1L])
  }
  call <- quote(npregbw(xdat = stop("not data yet"), regtype = "lp",
                        degree.select = "coordinate", search.engine = "nomad"))
  expect_true(.np_bw_call_uses_nomad_degree_search(call, caller))
  expect_false(.np_bw_call_uses_nomad_degree_search(quote(npregbw(regtype = "lc")), caller))
  expect_identical(hits, 0L)
  for (call in list(quote(user_selector(x)), quote(other::npregbw(x)),
                   as.call(list(function(x) x, quote(value)))))
    expect_identical(.np_bw_package_selector_call(call), call)
})

test_that("public consumers preserve values under caller masking", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  ns <- environment(npudens)
  e <- new.env(parent = ns)
  set.seed(334)
  e$x <- data.frame(x = runif(30), u = factor(rep(0:1,15)))
  e$xc <- e$x[1L]; e$y <- sin(e$x$x) + rnorm(30, sd=.2)
  e$dy <- data.frame(y=e$y); e$z <- data.frame(z=runif(30))
  e$d <- data.frame(e$x, y=e$y, z=e$z$z)
  cases <- list(
    reg=quote(npreg(txdat=xc,tydat=y,bws=.4)),
    density=quote(npudens(tdat=x,bws=c(.4,.2))),
    distribution=quote(npudist(tdat=xc,bws=.4)),
    cdensity=quote(npcdens(txdat=x,tydat=dy,bws=c(.4,.2,.4))),
    cdistribution=quote(npcdist(txdat=x,tydat=dy,bws=c(.4,.2,.4))),
    index=quote(npindex(txdat=xc,tydat=y,bws=c(1,.4),se=FALSE)),
    plreg=quote(npplreg(txdat=x,tydat=y,tzdat=z,bws=matrix(.4,3,1))),
    scoef=quote(npscoef(txdat=xc,tydat=y,tzdat=z,bws=.4)),
    qreg=quote(npqreg(txdat=x,tydat=dy,bws=c(.4,.2,.4))),
    conmode=quote(npconmode(txdat=x,tydat=data.frame(y=x$u),bws=c(.4,.2,.2))),
    sigtest=quote(npsigtest(xdat=x,ydat=y,bws=c(.4,.2),B=9,random.seed=3)),
    lsqreg=quote(nplsqreg(txdat=x,tydat=y,bws=c(.4,.2),bandwidth.compute=FALSE)),
    formula=quote(npudens(~x+u,data=d,bws=c(.4,.2))))
  fields <- function(out) list(fitted=if(!is.null(out$P))out$P else fitted(out),
    statistic=out$In, nobs=out$nobs)
  expected <- lapply(cases,function(call){set.seed(312);fields(eval(call,e))})
  for(name in c("npregbw","npplregbw","npudensbw","npudistbw","npcdensbw",
               "npcdistbw","npindexbw","npscoefbw"))
    makeActiveBinding(name,function(value)stop("masked selector was evaluated"),e)
  for(id in names(cases)) {
    set.seed(312);actual<-fields(eval(cases[[id]],e))
    expect_equal(actual,expected[[id]],tolerance=1e-12,info=id)
  }
})
