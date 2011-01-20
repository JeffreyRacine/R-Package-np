.onAttach <- function (lib, pkg) {
  cat("Nonparametric Kernel Methods for Mixed Datatypes (version 0.40-4)\n");
  if(is.null(options('np.messages')$np.messages))
    options(np.messages = TRUE)

  if(is.null(options('np.tree')$np.tree))
    options(np.tree = FALSE)

}
.Last.lib <- function (lpath){
  library.dynam.unload("np", libpath=lpath) 
  # cat("np unloaded\n")
}
