.onAttach <- function (lib, pkg) {
	packageStartupMessage("Nonparametric Kernel Methods for Mixed Datatypes (version 0.40-13)", domain = NULL,  appendLF = TRUE)

  if(is.null(options('np.messages')$np.messages))
    options(np.messages = TRUE)

  if(is.null(options('np.tree')$np.tree))
    options(np.tree = FALSE)

}
.Last.lib <- function (lpath){
  library.dynam.unload("np", libpath=lpath) 
}
