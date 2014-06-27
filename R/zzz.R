.onAttach <- function (lib, pkg) {
	packageStartupMessage("Nonparametric Kernel Methods for Mixed Datatypes (version 0.60-2)\n[vignette(\"np_faq\",package=\"np\") provides answers to frequently asked questions]", domain = NULL,  appendLF = TRUE)

  if(is.null(options('np.messages')$np.messages))
    options(np.messages = TRUE)

  if(is.null(options('np.tree')$np.tree))
    options(np.tree = FALSE)

}
.onUnload <- function (lpath){
  library.dynam.unload("np", libpath=lpath) 
}
