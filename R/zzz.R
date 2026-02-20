.onAttach <- function (lib, pkg) {
		packageStartupMessage("Parallel Nonparametric Kernel Methods for Mixed Datatypes (version 0.70-0) + Rmpi 0.7-3.3\n[vignette(\"np_faq\",package=\"npRmpi\") provides answers to frequently asked questions]\n[vignette(\"npRmpi\",package=\"npRmpi\") an overview]\n[vignette(\"entropy_np\",package=\"npRmpi\") an overview of entropy-based methods]", domain = NULL,  appendLF = TRUE)
    if (isTRUE(getOption("npRmpi.conflicts.warn", TRUE)) &&
        ("package:np" %in% search()) &&
        !isTRUE(getOption("npRmpi.conflicts.warned", FALSE))) {
      packageStartupMessage("note: both 'npRmpi' and 'np' are attached; prefer explicit npRmpi:: calls to avoid masking ambiguity")
      options(npRmpi.conflicts.warned = TRUE)
    }
}

.onUnload <- function (lpath){
  try(.C("np_release_static_buffers", as.integer(0), PACKAGE = "npRmpi"), silent = TRUE)
  mpi.finalize()
  library.dynam.unload("npRmpi", libpath=lpath) 
}

.onLoad <- function (lib, pkg) {
  library.dynam("npRmpi", pkg, lib)
  if (!TRUE)
    stop("Fail to load npRmpi dynamic library.")
  if (!is.loaded("mpi_initialize"))
    stop("Probably npRmpi has been detached. Please quit R.")

  # On macOS with MPICH, repeated spawn/close cycles have been observed to
  # destabilize subsequent MPI collectives. Default to reusing spawned slaves
  # (keep them alive) unless the user explicitly disables it.
  if (is.null(getOption("npRmpi.reuse.slaves")) &&
      identical(Sys.info()[["sysname"]], "Darwin") &&
      !nzchar(Sys.getenv("NP_RMPI_NO_REUSE_SLAVES"))) {
    options(npRmpi.reuse.slaves = TRUE)
  }

  if (nzchar(Sys.getenv("NP_RMPI_SKIP_INIT"))) {
    if(is.null(options('np.messages')$np.messages))
      options(np.messages = TRUE)

    if(is.null(options('np.tree')$np.tree))
      options(np.tree = FALSE)
    if(is.null(options('np.largeh.rel.tol')$np.largeh.rel.tol))
      options(np.largeh.rel.tol = 1e-3)
    if(is.null(options('np.disc.upper.rel.tol')$np.disc.upper.rel.tol))
      options(np.disc.upper.rel.tol = 1e-2)
    if(is.null(options('np.groupcv.fast')$np.groupcv.fast))
      options(np.groupcv.fast = TRUE)
    if (is.null(getOption("npRmpi.autodispatch")))
      options(npRmpi.autodispatch = FALSE)
    if (is.null(getOption("npRmpi.autodispatch.strict")))
      options(npRmpi.autodispatch.strict = TRUE)
    if (is.null(getOption("npRmpi.autodispatch.disable")))
      options(npRmpi.autodispatch.disable = FALSE)
    if (is.null(getOption("npRmpi.autodispatch.context")))
      options(npRmpi.autodispatch.context = FALSE)
    if (is.null(getOption("npRmpi.autodispatch.warned.nested")))
      options(npRmpi.autodispatch.warned.nested = FALSE)
    if (is.null(getOption("npRmpi.autodispatch.warned.bootstrap.plot")))
      options(npRmpi.autodispatch.warned.bootstrap.plot = FALSE)
    if (is.null(getOption("npRmpi.conflicts.warn")))
      options(npRmpi.conflicts.warn = TRUE)
    if (is.null(getOption("npRmpi.conflicts.warned")))
      options(npRmpi.conflicts.warned = FALSE)

    return(invisible())
  }

  if(.Call("mpidist",PACKAGE="npRmpi") == 2){
    auto.lamboot <- isTRUE(getOption("npRmpi.auto.lamboot", FALSE)) ||
      nzchar(Sys.getenv("NP_RMPI_AUTO_LAMBOOT"))
    if (auto.lamboot && (length(try(system("lamnodes",TRUE,ignore.stderr = TRUE))) == 0)){
	    system("lamboot -H",ignore.stderr = TRUE)
    }
  }
	
  if(!.Call("mpi_initialize",PACKAGE="npRmpi"))
    stop("Cannot start MPI_Init(). Exit")
  
  if(is.null(options('np.messages')$np.messages))
    options(np.messages = TRUE)

  if(is.null(options('np.tree')$np.tree))
    options(np.tree = FALSE)
  if(is.null(options('np.largeh.rel.tol')$np.largeh.rel.tol))
    options(np.largeh.rel.tol = 1e-3)
  if(is.null(options('np.disc.upper.rel.tol')$np.disc.upper.rel.tol))
    options(np.disc.upper.rel.tol = 1e-2)
  if(is.null(options('np.groupcv.fast')$np.groupcv.fast))
    options(np.groupcv.fast = TRUE)
  if (is.null(getOption("npRmpi.autodispatch")))
    options(npRmpi.autodispatch = FALSE)
  if (is.null(getOption("npRmpi.autodispatch.strict")))
    options(npRmpi.autodispatch.strict = TRUE)
  if (is.null(getOption("npRmpi.autodispatch.disable")))
    options(npRmpi.autodispatch.disable = FALSE)
  if (is.null(getOption("npRmpi.autodispatch.context")))
    options(npRmpi.autodispatch.context = FALSE)
  if (is.null(getOption("npRmpi.autodispatch.warned.nested")))
    options(npRmpi.autodispatch.warned.nested = FALSE)
  if (is.null(getOption("npRmpi.autodispatch.warned.bootstrap.plot")))
    options(npRmpi.autodispatch.warned.bootstrap.plot = FALSE)
  if (is.null(getOption("npRmpi.conflicts.warn")))
    options(npRmpi.conflicts.warn = TRUE)
  if (is.null(getOption("npRmpi.conflicts.warned")))
    options(npRmpi.conflicts.warned = FALSE)
}
