# Test harness only. Opt-in explicit masters avoid untracked singleton MPI
# launchers. Do not change package initialization, spawn or attach semantics.
npRmpi_test_launch_spec <- function(command, args=character(),
                                   launcher=Sys.getenv("NP_RMPI_TEST_MPIEXEC",""),
                                   inherited=names(Sys.getenv()),
                                   env=character()) {
  if (!nzchar(launcher))
    return(list(command=command,args=args))
  if (.Platform$OS.type != "unix" || !startsWith(launcher,"/") ||
      !file.exists(launcher) || file.access(launcher,1L) != 0L)
    stop("explicit test MPI launcher must be an executable absolute Unix path")
  bootstrap <- c("PMI_FD","PMI_PORT","PMI_RANK","PMI_SIZE","PMI_ID","PMI_HOSTNAME",
                 "PMI_APPNUM","PMI_KVS","PMI_KVSNAME",
                 "MPI_LOCALRANKID","MPI_LOCALNRANKS")
  bootstrap_family <- "^(PMI_|PMIX_|OMPI_COMM_WORLD_|OMPI_MCA_ess_|OMPI_MCA_orte_|MPI_LOCAL)"
  candidates <- inherited[grepl(bootstrap_family,inherited)]
  if(length(setdiff(candidates,bootstrap)))
    stop("unqualified inherited MPI bootstrap variables: ",
         paste(setdiff(candidates,bootstrap),collapse=", "))
  supplied <- sub("=.*$","",env)
  if(any(supplied %in% bootstrap) ||
     any(grepl(bootstrap_family,supplied)))
    stop("independent test child cannot override MPI bootstrap identity")
  if(identical(basename(command),"Rscript")) {
    args <- c("-n","1",shQuote(command),args)
    command <- launcher
  }
  remove <- intersect(inherited,bootstrap)
  unset <- unlist(lapply(remove,function(key) c("-u",key)),use.names=FALSE)
  list(command="/usr/bin/env",args=c(unset,shQuote(command),args))
}

npRmpi_test_system2 <- function(command, args=character(), ...,
                               env=character()) {
  launch <- npRmpi_test_launch_spec(command,args,env=env)
  system2(launch$command,launch$args,...,env=env)
}
