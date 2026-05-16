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
  if (is.null(getOption("np.categorical.compress")))
    options(np.categorical.compress = TRUE)
  if (is.null(getOption("np.largeh.rel.tol")))
    options(np.largeh.rel.tol = 1e-3)
  if (is.null(getOption("np.disc.upper.rel.tol")))
    options(np.disc.upper.rel.tol = 1e-2)
}

npTreeOrCategoricalCompress <- function(ncon = 0L, ncat = 0L) {
  if (isTRUE(getOption("np.tree")))
    return(TRUE)

  if (!isTRUE(getOption("np.categorical.compress", TRUE)))
    return(FALSE)

  ncon <- if (is.null(ncon) || !length(ncon)) 0L else sum(as.integer(ncon), na.rm = TRUE)
  ncat <- if (is.null(ncat) || !length(ncat)) 0L else sum(as.integer(ncat), na.rm = TRUE)

  isTRUE(ncon == 0L) && isTRUE(ncat > 0L)
}

npDoTreeOrCategoricalCompress <- function(ncon = 0L, ncat = 0L) {
  if (npTreeOrCategoricalCompress(ncon = ncon, ncat = ncat)) {
    DO_TREE_YES
  } else {
    DO_TREE_NO
  }
}

.onUnload <- function (lpath){
  tryCatch(.Call("C_np_release_static_buffers", PACKAGE = "np"), error = function(e) NULL)
  library.dynam.unload("np", libpath=lpath) 
}
