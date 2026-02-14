.onAttach <- function (lib, pkg) {
	packageStartupMessage("Nonparametric Kernel Methods for Mixed Datatypes (version 0.60-21)\n[vignette(\"np_faq\",package=\"np\") provides answers to frequently asked questions]\n[vignette(\"np\",package=\"np\") an overview]\n[vignette(\"entropy_np\",package=\"np\") an overview of entropy-based methods]", domain = NULL,  appendLF = TRUE)
}

.onLoad <- function (lib, pkg) {
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
}

.onUnload <- function (lpath){
  try(.C("np_release_static_buffers", as.integer(0), PACKAGE = "np"), silent = TRUE)
  library.dynam.unload("np", libpath=lpath) 
}
