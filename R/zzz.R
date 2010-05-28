.onAttach <- function (lib, pkg) {
  cat("Nonparametric Kernel Methods for Mixed Datatypes (version 0.40-1)\n");
  if(is.null(options('np.messages')$np.messages))
    options(np.messages = TRUE)

}
.Last.lib <- function (lpath){
  library.dynam.unload("np", libpath=lpath) 
  # cat("np unloaded\n")
}
