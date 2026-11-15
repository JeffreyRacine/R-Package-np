test_that("quantile progress inspection never executes materialized options", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI session unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages=FALSE)
  on.exit(options(old),add=TRUE)
  d <- data.frame(x=seq(-1,1,length.out=30L),y=sin(seq_len(30L)/4))
  flag <- new.env(parent=emptyenv())
  for(route in c("positional","named"))
    for(control in c("regtype","degree.select","search.engine")) {
      flag$n <- 0L
      bad <- switch(control,
        regtype=quote({flag$n <- flag$n+1L; "lc"}),
        degree.select=quote({flag$n <- flag$n+1L; "coordinate"}),
        search.engine=quote({flag$n <- flag$n+1L; "nomad"}))
      args <- c(if(route=="positional") list(bws=y~x) else list(formula=y~x),
        list(data=d,nomad=FALSE,regtype="lp",degree.select="coordinate",
             search.engine="nomad",bandwidth.compute=FALSE))
      args[[control]] <- bad
      set.seed(8341); seed <- .Random.seed
      expect_error(do.call(npqreg,args,quote=TRUE),"executable option",fixed=TRUE)
      expect_identical(flag$n,0L)
      expect_identical(.Random.seed,seed)
    }
})

test_that("literal option inspection preserves degree-search progress choices", {
  inspector <- getFromNamespace(".np_bw_call_uses_nomad_degree_search","npRmpi")
  seen <- NULL
  local_mocked_bindings(
    .np_bw_call_uses_nomad_degree_search=function(call_obj,caller_env) {
      seen <<- inspector(call_obj,caller_env);seen
    },
    .np_formula_dispatch_call=function(...) stop("inspection completed"),
    .np_progress_select_bandwidth_enhanced=function(label,expr) force(expr),
    .package="npRmpi")
  method <- getS3method("npqreg","default")
  d <- data.frame(x=seq_len(6L),y=seq_len(6L))
  omit <- structure(list(),class="omit")
  controls <- list(regtype=list("lc","ll","lp",omit,NULL),
    select=list("manual","coordinate","exhaustive",omit,NULL),
    engine=list("nomad+powell","cell","nomad",omit,NULL))
  resolved <- function(x,default) if(inherits(x,"omit")||is.null(x)) default else x
  for(route in c("positional","named")) for(nomad in list(FALSE,TRUE,"auto"))
    for(regtype in controls$regtype) for(select in controls$select)
      for(engine in controls$engine) {
        args <- c(if(route=="positional")list(bws=y~x)else list(formula=y~x),
          list(data=d,nomad=nomad))
        vals <- list(regtype=regtype,degree.select=select,search.engine=engine)
        for(name in names(vals)) if(!inherits(vals[[name]],"omit"))
          args[name] <- vals[name]
        seen <- NULL
        expect_error(do.call(method,args),"inspection completed",fixed=TRUE)
        expected <- !identical(nomad,FALSE) ||
          (resolved(regtype,"lc")=="lp" && resolved(select,"manual")!="manual" &&
           resolved(engine,"nomad+powell") %in% c("nomad","nomad+powell"))
        expect_identical(seen,expected)
      }
})
