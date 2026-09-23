test_that("package call ownership preserves explicit heads and custom definitions", {
  ns <- asNamespace("npRmpi")
  bind <- get(".npRmpi_autodispatch_bind_call_owner",ns)
  f <- get("npcmstest",ns)
  mc <- quote(local_alias(xdat=side_effect()))
  out <- bind(mc,"npcmstest",mc,f)
  expect_true(attr(out,".npRmpi.bound.call.owner",exact=TRUE))
  attr(out,".npRmpi.bound.call.owner") <- NULL
  expect_identical(out,quote(npcmstest(xdat=side_effect())))
  for(head in list(quote(factory()),quote(holder$fun),f,quote(get("fun")))) {
    mc[[1L]] <- head
    expect_identical(bind(mc,"npcmstest",mc,f)[[1L]],as.name("npcmstest"))
  }
  rewritten <- mc;rewritten[[1L]] <- get(".npRmpi_lease_context_symbol",ns)()
  expect_identical(bind(rewritten,"npcmstest",mc,f),rewritten)
  expect_identical(bind(mc,"npcmstest",mc,function(...)NULL),mc)
  expect_identical(bind(mc,NULL,mc,f),mc)
  canonical <- quote(npcmstest(xdat=side_effect()))
  expect_identical(bind(canonical,"npcmstest",canonical,function(...)NULL),canonical)
  method <- get("npregbw.default",ns)
  generic <- quote(npregbw(xdat=side_effect()))
  bound <- bind(generic,"npregbw.default",quote(npregbw.default(xdat=side_effect())),method)
  expect_identical(bound[[1L]],as.name("npregbw.default"))
  expect_true(attr(bound,".npRmpi.bound.call.owner",exact=TRUE))
})

test_that("MPI package aliases evaluate function heads and data exactly once", {
  skip_on_cran()
  if(!spawn_mpi_slaves(1L))skip("MPI slaves unavailable")
  on.exit(close_mpi_slaves(),add=TRUE)
  old <- options(np.messages=FALSE);on.exit(options(old),add=TRUE)
  d <- data.frame(x=seq(-1,1,length.out=31))
  d$y <- d$x+sin(seq_len(31))*.2
  model <- lm(y~x,data=d,x=TRUE,y=TRUE)
  args <- list(xdat=d["x"],ydat=d$y,model=model,bws=.4,
    bandwidth.compute=FALSE,distribution="asymptotic")
  ref <- do.call(npcmstest,args)
  payload <- function(z) unlist(z[c("Tn","P","Jn","Omega")])
  expect_true(length(payload(ref))>0)
  for(shape in c("alias","list","get","factory","do.call","wrapper")) {
    hits <- 0L;data.hits <- 0L
    fun <- npcmstest;holder <- list(fun=fun)
    factory <- function() {hits<<-hits+1L;fun}
    data <- function() {data.hits<<-data.hits+1L;d["x"]}
    expr <- switch(shape,
      alias=quote(fun(xdat=data(),ydat=d$y,model=model,bws=.4,
        bandwidth.compute=FALSE,distribution="asymptotic")),
      list=quote(holder$fun(xdat=data(),ydat=d$y,model=model,bws=.4,
        bandwidth.compute=FALSE,distribution="asymptotic")),
      get=quote(get("fun")(xdat=data(),ydat=d$y,model=model,bws=.4,
        bandwidth.compute=FALSE,distribution="asymptotic")),
      factory=quote(factory()(xdat=data(),ydat=d$y,model=model,bws=.4,
        bandwidth.compute=FALSE,distribution="asymptotic")),
      do.call=quote(do.call(fun,c(args[setdiff(names(args),"xdat")],list(xdat=data())))),
      wrapper=quote((function(...) fun(...))(xdat=data(),ydat=d$y,model=model,bws=.4,
        bandwidth.compute=FALSE,distribution="asymptotic")))
    z <- eval(expr)
    expect_equal(payload(z),payload(ref),tolerance=0,info=shape)
    expect_identical(data.hits,1L,info=shape)
    expect_identical(hits,as.integer(shape=="factory"),info=shape)
  }
  # A shadowed local canonical symbol cannot replace the selected alias.
  z <- local({npcmstest<-function(...)stop("wrong function");fun<-npRmpi::npcmstest
    fun(xdat=d["x"],ydat=d$y,model=model,bws=.4,
      bandwidth.compute=FALSE,distribution="asymptotic")})
  expect_equal(payload(z),payload(ref),tolerance=0)
  mpi.bcast.cmd(assign("npcmstest",function(...)stop("wrong global function"),
                      envir=.GlobalEnv),caller.execute=TRUE)
  tryCatch({
    z <- fun(xdat=d["x"],ydat=d$y,model=model,bws=.4,
      bandwidth.compute=FALSE,distribution="asymptotic")
    expect_equal(payload(z),payload(ref),tolerance=0)
  },finally=mpi.bcast.cmd(rm("npcmstest",envir=.GlobalEnv),caller.execute=TRUE))
  bwfun <- npregbw
  b <- bwfun(y~x,data=d,bws=.4,bandwidth.compute=FALSE)
  expect_true(is.symbol(b$call[[1L]]))
  expect_null(attr(b$call,".npRmpi.bound.call.owner",exact=TRUE))
  mpi.bcast.cmd(assign("npregbw",function(...)stop("wrong global constructor"),
                      envir=.GlobalEnv),caller.execute=TRUE)
  tryCatch({
    masked <- bwfun(y~x,data=d,bws=.4,bandwidth.compute=FALSE)
    expect_identical(masked$bw,b$bw)
  },finally=mpi.bcast.cmd(rm("npregbw",envir=.GlobalEnv),caller.execute=TRUE))
  fitfun <- npreg
  expect_equal(fitted(fitfun(b)),fitted(npreg(b)),tolerance=0)
  expect_equal(fitted(fitfun(y~x,data=d,bws=.4)),fitted(npreg(b)),tolerance=0)
  expect_error(fun(xdat=d["x"],ydat=d$y,model=model,bws=.4,
    bandwidth.compute=FALSE,distribution="bad"),"arg")
})
