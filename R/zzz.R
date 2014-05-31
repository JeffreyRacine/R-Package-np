.onAttach <- function (lib, pkg) {
	packageStartupMessage("Nonparametric Kernel Methods for Mixed Datatypes (version 0.60-0) + Rmpi 0.6-5\n[vignette(\"np_faq\",package=\"npRmpi\") provides answers to frequently asked questions]", domain = NULL,  appendLF = TRUE)

  if(is.null(options('np.messages')$np.messages))
    options(np.messages = TRUE)

  if(is.null(options('np.tree')$np.tree))
    options(np.tree = FALSE)
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
}
