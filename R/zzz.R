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
  if (is.null(getOption("np.powell.cache")))
    options(np.powell.cache = .np_powell_cache_default())
  if (is.null(getOption("np.largeh")))
    options(np.largeh = TRUE)
  if (is.null(getOption("np.largelambda")))
    options(np.largelambda = TRUE)
  if (is.null(getOption("np.largeh.rel.tol")))
    options(np.largeh.rel.tol = 1e-3)
  if (is.null(getOption("np.disc.upper.rel.tol")))
    options(np.disc.upper.rel.tol = 1e-2)
}

.np_powell_cache_default <- function() {
  flag <- Sys.getenv("NP_NN_POWELL_CACHE_INSTRUMENT", unset = "")
  if (!nzchar(flag))
    return(TRUE)
  !(flag %in% c("0", "false", "FALSE", "off", "OFF", "no", "NO"))
}

npCountVars <- function(x) {
  if (is.null(x) || !length(x))
    return(0L)
  sum(as.integer(x), na.rm = TRUE)
}

npUseContinuousTree <- function(ncon = 0L) {
  isTRUE(getOption("np.tree")) && isTRUE(npCountVars(ncon) > 0L)
}

npUseCategoricalCompress <- function(ncon = 0L, ncat = 0L) {
  isTRUE(getOption("np.categorical.compress", TRUE)) &&
    isTRUE(npCountVars(ncon) == 0L) &&
    isTRUE(npCountVars(ncat) > 0L)
}

npUseKernelAccelerationFlag <- function(ncon = 0L, ncat = 0L) {
  npUseContinuousTree(ncon = ncon) ||
    npUseCategoricalCompress(ncon = ncon, ncat = ncat)
}

npDoTreeOrCategoricalCompress <- function(ncon = 0L, ncat = 0L) {
  if (npUseKernelAccelerationFlag(ncon = ncon, ncat = ncat)) {
    DO_TREE_YES
  } else {
    DO_TREE_NO
  }
}

.onUnload <- function (lpath){
  tryCatch(.Call("C_np_release_static_buffers", PACKAGE = "np"), error = function(e) NULL)
  library.dynam.unload("np", libpath=lpath) 
}
