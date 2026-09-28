test_that("conditional CDF handoff retains executed Powell cache reporting", {
  old <- options(np.messages=FALSE,np.tree=FALSE,np.objective.cache=TRUE)
  on.exit(options(old),add=TRUE)
  pkg <- getNamespaceName(environment(npcdistbw)); ns <- asNamespace(pkg)
  capture <- new.env(parent=emptyenv()); capture$hot <- list()
  trace(".npcdistbw_run_fixed_degree",where=ns,print=FALSE,
    exit=substitute({
      value <- returnValue()
      assign("hot",c(get("hot",CAPTURE),list(list(cache=value$nn.cache,n=value$num.feval))),envir=CAPTURE)
    },list(CAPTURE=capture)))
  trace(".npcdistbw_run_fixed_degree_mads",where=ns,print=FALSE,
    exit=substitute(assign("native",mads.num.feval.total,envir=CAPTURE),list(CAPTURE=capture)))
  trace(".npcdistbw_nomad_search",where=ns,print=FALSE,
    exit=substitute(assign("native",nomad.num.feval.total,envir=CAPTURE),list(CAPTURE=capture)))
  on.exit({
    untrace(".npcdistbw_run_fixed_degree",where=ns)
    untrace(".npcdistbw_run_fixed_degree_mads",where=ns)
    untrace(".npcdistbw_nomad_search",where=ns)
  },add=TRUE)
  x <- data.frame(x=seq(-1,1,length.out=24L)); y <- data.frame(y=sin(2*x$x)+seq_len(24)/100)
  for(case in c("gnn","ann","disabled","restarts","degree","mads-only")) {
    capture$hot <- list(); capture$native <- NA_real_; options(np.objective.cache=case!="disabled")
    args <- list(xdat=x,ydat=y,regtype="lc",bwtype=if(case=="gnn") "generalized_nn" else "adaptive_nn",
      bwsolver=if(case=="mads-only") "mads" else "mads+powell",nmulti=if(case=="restarts") 2L else 1L,
      itmax=120L,powell.remin=FALSE,nomad.opts=list(MAX_BB_EVAL=40L))
    if(case=="degree") {
      args$bwsolver <- NULL; args$search.engine <- "nomad+powell"
      args$regtype <- "lp"; args$nomad <- TRUE; args$degree.min <- 0L; args$degree.max <- 1L
    }
    b <- do.call(npcdistbw,args)
    native <- as.numeric(capture$native)
    if(case=="mads-only") {
      expect_length(capture$hot,0L); expect_identical(as.numeric(b$num.feval),native)
    } else {
      expect_length(capture$hot,1L)
      hot <- capture$hot[[1L]]
      expect_identical(b$nn.cache,hot$cache,info=case)
      expect_identical(as.numeric(b$num.feval),native+as.numeric(hot$n),info=case)
      if(case=="disabled") expect_true(all(hot$cache==0)) else {
        expect_identical(as.numeric(sum(hot$cache[c("raw.evals","hits")])),as.numeric(hot$n))
        expect_gt(hot$cache[["hits"]],0)
      }
    }
    raw <- getFromNamespace(".npcdistbw_eval_only",pkg)(xdat=x,ydat=y,bws=b,invalid.penalty="dbmax")$objective
    expect_identical(as.numeric(raw),as.numeric(b$fval),info=case)
  }
})
