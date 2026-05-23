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
    options(np.tree = "auto")
  if (is.null(getOption("np.categorical.compress")))
    options(np.categorical.compress = TRUE)
  if (is.null(getOption("np.powell.cache")))
    options(np.powell.cache = .np_powell_cache_default())
  if (is.null(getOption("np.largeh")))
    options(np.largeh = TRUE)
  if (is.null(getOption("np.largelambda")))
    options(np.largelambda = TRUE)
  if (is.null(getOption("np.largenn")))
    options(np.largenn = FALSE)
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

npTreeMode <- function(value = getOption("np.tree", "auto")) {
  if (isTRUE(value))
    return("on")
  if (identical(value, FALSE))
    return("off")
  if (is.character(value) && length(value) == 1L &&
      identical(unname(tolower(value)), "auto"))
    return("auto")
  stop("option 'np.tree' must be TRUE, FALSE, or \"auto\"", call. = FALSE)
}

npTreeContinuousKernelTypes <- function(bws = NULL, ckertype = NULL) {
  if (!is.null(ckertype))
    return(as.character(ckertype))
  if (is.null(bws))
    return(character(0L))

  out <- character(0L)
  if (!is.null(bws$yncon) || !is.null(bws$xncon)) {
    if (npCountVars(bws$yncon) > 0L && !is.null(bws$cykertype))
      out <- c(out, as.character(bws$cykertype))
    if (npCountVars(bws$xncon) > 0L && !is.null(bws$cxkertype))
      out <- c(out, as.character(bws$cxkertype))
    return(out)
  }

  if (npCountVars(bws$ncon) > 0L && !is.null(bws$ckertype))
    return(as.character(bws$ckertype))
  character(0L)
}

npTreeAutoKernelEligible <- function(bws = NULL, ckertype = NULL) {
  kernels <- npTreeContinuousKernelTypes(bws = bws, ckertype = ckertype)
  length(kernels) > 0L &&
    all(tolower(kernels) %in% c("epanechnikov", "uniform"))
}

npUseContinuousTree <- function(ncon = 0L, bws = NULL, ckertype = NULL) {
  if (!isTRUE(npCountVars(ncon) > 0L))
    return(FALSE)

  mode <- npTreeMode()
  if (identical(mode, "on"))
    return(TRUE)
  if (identical(mode, "off"))
    return(FALSE)

  npTreeAutoKernelEligible(bws = bws, ckertype = ckertype)
}

npUseCategoricalCompress <- function(ncon = 0L, ncat = 0L) {
  isTRUE(getOption("np.categorical.compress", TRUE)) &&
    isTRUE(npCountVars(ncon) == 0L) &&
    isTRUE(npCountVars(ncat) > 0L)
}

npUseKernelAccelerationFlag <- function(ncon = 0L, ncat = 0L,
                                        bws = NULL, ckertype = NULL) {
  npUseContinuousTree(ncon = ncon, bws = bws, ckertype = ckertype) ||
    npUseCategoricalCompress(ncon = ncon, ncat = ncat)
}

npDoTreeOrCategoricalCompress <- function(ncon = 0L, ncat = 0L,
                                          bws = NULL, ckertype = NULL) {
  if (npUseKernelAccelerationFlag(ncon = ncon, ncat = ncat,
                                  bws = bws, ckertype = ckertype)) {
    DO_TREE_YES
  } else {
    DO_TREE_NO
  }
}

.onUnload <- function (lpath){
  tryCatch(.Call("C_np_release_static_buffers", PACKAGE = "np"), error = function(e) NULL)
  library.dynam.unload("np", libpath=lpath) 
}
