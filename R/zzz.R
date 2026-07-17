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
  if (is.null(getOption("np.objective.cache")))
    options(np.objective.cache = TRUE)
  if (is.null(getOption("np.largeh")))
    options(np.largeh = TRUE)
  if (is.null(getOption("np.largelambda")))
    options(np.largelambda = TRUE)
  if (is.null(getOption("np.extendednn")))
    options(np.extendednn = TRUE)
  if (is.null(getOption("np.macMseries.accelerate")))
    options(np.macMseries.accelerate = "auto")
  if (is.null(getOption("np.nomad.degree.start.policy")))
    options(np.nomad.degree.start.policy = "low_first_full_random")
  if (is.null(getOption("np.largeh.rel.tol")))
    options(np.largeh.rel.tol = 1e-3)
  if (is.null(getOption("np.disc.upper.rel.tol")))
    options(np.disc.upper.rel.tol = 1e-2)
}

npCountVars <- function(x) {
  if (is.null(x) || !length(x))
    return(0L)
  sum(as.integer(x), na.rm = TRUE)
}

npObjectiveCacheEnabled <- function(value = getOption("np.objective.cache", TRUE)) {
  if (isTRUE(value))
    return(TRUE)
  if (identical(value, FALSE))
    return(FALSE)
  stop("option 'np.objective.cache' must be TRUE or FALSE", call. = FALSE)
}

npLogicalOption <- function(name, default = TRUE) {
  value <- getOption(name, default)
  if (isTRUE(value))
    return(TRUE)
  if (identical(value, FALSE))
    return(FALSE)
  stop(sprintf("option '%s' must be TRUE or FALSE", name), call. = FALSE)
}

npFastpathTolerance <- function(name, default, lower, upper) {
  value <- getOption(name, default)
  if (is.numeric(value) && length(value) == 1L && !is.na(value) &&
      is.finite(value) && value > lower && value < upper)
    return(as.double(value))
  stop(
    sprintf("option '%s' must be a finite numeric scalar in (%s, %s)",
            name, format(lower), format(upper)),
    call. = FALSE
  )
}

npLargehRelTol <- function() {
  npFastpathTolerance("np.largeh.rel.tol", 1e-3, 0, 0.1)
}

npDiscUpperRelTol <- function() {
  npFastpathTolerance("np.disc.upper.rel.tol", 1e-2, 0, 0.5)
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
