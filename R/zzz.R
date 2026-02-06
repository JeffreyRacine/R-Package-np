.onAttach <- function (lib, pkg) {
	packageStartupMessage("Parallel Nonparametric Kernel Methods for Mixed Datatypes (version 0.60-20) + Rmpi 0.7-3.3\n[vignette(\"np_faq\",package=\"npRmpi\") provides answers to frequently asked questions]\n[vignette(\"npRmpi\",package=\"npRmpi\") an overview]\n[vignette(\"entropy_np\",package=\"npRmpi\") an overview of entropy-based methods]", domain = NULL,  appendLF = TRUE)
}

.onUnload <- function (lpath){
  mpi.finalize()
  library.dynam.unload("npRmpi", libpath=lpath) 
}

.onLoad <- function (lib, pkg) {
  library.dynam("npRmpi", pkg, lib)
  if (!TRUE)
    stop("Fail to load npRmpi dynamic library.")
  if (!is.loaded("mpi_initialize"))
    stop("Probably npRmpi has been detached. Please quit R.")

  if (nzchar(Sys.getenv("NP_RMPI_SKIP_INIT"))) {
    if(is.null(options('np.messages')$np.messages))
      options(np.messages = TRUE)

    if(is.null(options('np.tree')$np.tree))
      options(np.tree = FALSE)

    return(invisible())
  }

  if(.Call("mpidist",PACKAGE="npRmpi") == 2){
    if (length(try(system("lamnodes",TRUE,ignore.stderr = TRUE))) == 0){
	    system("lamboot -H",ignore.stderr = TRUE)
    }
  }
	
  if(!.Call("mpi_initialize",PACKAGE="npRmpi"))
    stop("Cannot start MPI_Init(). Exit")
  
  if (exists(".Random.seed") && 
      round(.Random.seed[1]-5,-1) == .Random.seed[1]-5) {
    rm(.Random.seed, envir=.GlobalEnv)
  }

  if(is.null(options('np.messages')$np.messages))
    options(np.messages = TRUE)

  if(is.null(options('np.tree')$np.tree))
    options(np.tree = FALSE)
}
