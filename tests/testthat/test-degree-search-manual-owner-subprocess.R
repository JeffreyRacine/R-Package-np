test_that("manual search uses the same collective owner on every rank", {
  skip_if(Sys.getenv("NP_RMPI_RUN_SEARCH_OWNER_TESTS") != "true",
          "explicit MPI search-owner qualification")
  env <- npRmpi_subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable")
  for (workers in c(1L, 2L)) {
    res <- npRmpi_run_rscript_subprocess(lines = c(
      "library(npRmpi)",
      sprintf("npRmpi.init(nslaves = %dL, quiet = TRUE)", workers),
      "mpi.bcast.cmd({options(np.messages=FALSE); assign('.r40_seen', logical(), .GlobalEnv); trace('.npRmpi_degree_search_collective', where=asNamespace('npRmpi'), exit=quote(.GlobalEnv$.r40_seen <- c(.GlobalEnv$.r40_seen, returnValue())), print=FALSE); set.seed(42); x<-seq(-1,1,length.out=30); y<-sin(x)}, caller.execute=TRUE)",
      "for (family in c('npreg', 'npcdens', 'npcdist')) {",
      "  mc <- substitute(FUN(y~x, nomad=TRUE, degree.max=1L, nmulti=1L), list(FUN=as.name(family)))",
      "  fit <- get('.npRmpi_bcast_cmd_expr',asNamespace('npRmpi'))(mc, caller.execute=TRUE)",
      "  stopifnot(all(is.finite(fitted(fit))))",
      "}",
      "owners <- mpi.bcast.cmd(get('mpi.gather.Robj', asNamespace('npRmpi'))(.GlobalEnv$.r40_seen, simplify=FALSE), caller.execute=TRUE)",
      "stopifnot(all(vapply(owners, function(x) length(x)>=3L && all(x), TRUE)))",
      "npRmpi.quit(force=TRUE); mpi.quit()"
    ), timeout=120L, env=env)
    expect_equal(res$status, 0L, info=paste(res$output, collapse="\n"))
  }
})
