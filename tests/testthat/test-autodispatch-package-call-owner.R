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

test_that("dots expansion retains only the existing call-owner marker", {
 f<-get(".npRmpi_autodispatch_expand_dots_call",asNamespace("npRmpi"))
 for(marker in list(NULL,FALSE,TRUE)) for(dots in list(pairlist(),pairlist(a=1,2),list(a=1,2),quote(list(a=1,2)),quote(pairlist(a=1,2)))) {
   x<-as.call(c(list(as.name("selected")),list(x=quote(side_effect()),"..."=dots,z=3)))
   attr(x,".npRmpi.bound.call.owner")<-marker
   y<-f(x);expect_identical(attr(y,".npRmpi.bound.call.owner",exact=TRUE),marker)
   attr(y,".npRmpi.bound.call.owner")<-NULL
   expected<-if(length(dots)==0L)quote(selected(x=side_effect(),z=3)) else quote(selected(x=side_effect(),a=1,2,z=3))
   expect_identical(y,expected)
 }
 for(x in list(quote(selected()),quote(selected(x=side_effect())),NULL)) {
   if(!is.null(x))attr(x,".npRmpi.bound.call.owner")<-TRUE
   expect_identical(f(x),x)
 }
})

test_that("fresh namespace-only aliases preserve bandwidth call ownership", {
  skip_on_cran()
  lines <- c(
    "ns<-loadNamespace(\"npRmpi\")",
    "options(np.messages=FALSE)",
    "stopifnot(!\"package:npRmpi\"%in%search())",
    "ns$npRmpi.init(nslaves=as.integer(Sys.getenv(\"N5_SLAVES\",\"1\")),quiet=TRUE)",
    "on.exit(ns$npRmpi.quit(force=TRUE),add=TRUE)",
    "set.seed(811);d<-data.frame(y=rnorm(80),a=rnorm(80),b=runif(80))",
    "cases<-list(sc=list(fun=\"npscoefbw\",expr=quote(f(y~a|b,data=d,nmulti=1))),",
    " pl=list(fun=\"npplregbw\",expr=quote(f(y~a|b,data=d,nmulti=1))),",
    " ix=list(fun=\"npindexbw\",expr=quote(f(y~a+b,data=d,nmulti=1))))",
    "if(!identical(Sys.getenv(\"N5_FULL\"),\"true\"))cases<-cases[1]",
    "values<-list()",
    "for(n in names(cases)){",
    " f<-get(cases[[n]]$fun,ns);cat(\"BEGIN\",n,\"\\n\");flush(stdout())",
    " set.seed(811);v<-eval(cases[[n]]$expr)",
    " values[[n]]<-v;stopifnot(is.list(v));cat(\"DONE\",n,\"\\n\");flush(stdout())",
    "}",
    "library(npRmpi)",
    "for(n in names(cases)){",
    " f<-get(cases[[n]]$fun,ns);set.seed(811);v<-eval(cases[[n]]$expr)",
    " for(k in intersect(c(\"bw\",\"xbw\",\"ybw\",\"xnames\",\"ynames\",\"znames\"),names(v))) {",
    "   a<-v[[k]];b<-values[[n]][[k]]",
    "   if(k==\"bw\" && is.list(a)){a<-lapply(a,function(x)x$bw);b<-lapply(b,function(x)x$bw)}",
    "   if(!identical(a,b)){cat(\"MISMATCH\",n,k,\"\\n\");print(all.equal(a,b));stop(\"public payload differs\")}",
    " }",
    "}",
    "npRmpi.quit(force=TRUE)",
    "cat(\"N5_OWNER_NAMESPACE_COMPLETE\\n\")"
  )
  out <- npRmpi_run_rscript_subprocess(lines, timeout=90L)
  expect_identical(out$status,0L)
  expect_true(any(out$output=="N5_OWNER_NAMESPACE_COMPLETE"))
})
