.onAttach <- function (lib, pkg) {
	packageStartupMessage(
    sprintf(
      paste(
        "np %s",
        "Examples and guides at https://jeffreyracine.github.io/gallery/",
        'See also vignette("np_getting_started", package = "np")',
        sep = "\n"
      ),
      utils::packageDescription(pkg, fields = "Version")
    ),
    domain = NULL,
    appendLF = TRUE
  )
}

.onLoad <- function (lib, pkg) {
  if (is.null(getOption("np.messages")))
    options(np.messages = TRUE)
  if (is.null(getOption("np.tree")))
    options(np.tree = FALSE)
  if (is.null(getOption("np.largeh.rel.tol")))
    options(np.largeh.rel.tol = 1e-3)
  if (is.null(getOption("np.disc.upper.rel.tol")))
    options(np.disc.upper.rel.tol = 1e-2)
}

.onUnload <- function (lpath){
  tryCatch(.Call("C_np_release_static_buffers", PACKAGE = "np"), error = function(e) NULL)
  library.dynam.unload("np", libpath=lpath) 
}
