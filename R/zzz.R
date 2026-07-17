.npRmpi_embedded_backend_version <- "0.7-3.3"

.npRmpi_s3_generic_namespace <- function(generic) {
  switch(
    generic,
    "plot" = "graphics",
    "predict" = "stats",
    "fitted" = "stats",
    "residuals" = "stats",
    "coef" = "stats",
    "vcov" = "stats",
    "quantile" = "stats",
    "print" = "base",
    "summary" = "base",
    "mode" = "base",
    "as.double" = "base",
    NULL
  )
}

.npRmpi_reassert_s3_methods <- function(generics = NULL) {
  ns <- asNamespace("npRmpi")
  s3 <- getNamespaceInfo(ns, "S3methods")

  if (is.null(s3) || !length(s3))
    return(invisible(TRUE))

  if (is.null(dim(s3)))
    s3 <- matrix(s3, ncol = 4L)

  if (is.null(generics)) {
    generics <- unique(as.character(s3[, 1L]))
  } else {
    generics <- unique(as.character(generics))
  }

  for (generic in generics) {
    generic.pkg <- .npRmpi_s3_generic_namespace(generic)
    if (is.null(generic.pkg))
      next

    generic.ns <- asNamespace(generic.pkg)
    rows <- which(as.character(s3[, 1L]) == generic)
    if (!length(rows))
      next

    for (row in rows) {
      cls <- as.character(s3[row, 2L])
      method.name <- as.character(s3[row, 3L])

      if (!nzchar(method.name) || !exists(method.name, envir = ns, inherits = FALSE))
        next

      registerS3method(
        generic,
        cls,
        get(method.name, envir = ns, inherits = FALSE),
        envir = generic.ns
      )
    }
  }

  invisible(TRUE)
}

.npRmpi_after_np_load <- function() {
  tryCatch(.npRmpi_reassert_s3_methods(), error = function(e) invisible(e))
  invisible(TRUE)
}

.npRmpi_after_np_attach <- function() {
  tryCatch(.npRmpi_after_np_load(), error = function(e) invisible(e))
  tryCatch(get(".npRmpi_warn_pkg_conflict_once",
               envir = asNamespace("npRmpi"),
               mode = "function",
               inherits = FALSE)(),
           error = function(e) invisible(e))
  invisible(TRUE)
}

.npRmpi_try_dynload <- function(lib, pkg) {
  dll <- paste0(pkg, .Platform$dynlib.ext)
  candidates <- unique(c(
    file.path(lib, pkg, "libs", dll),
    file.path(lib, "libs", dll),
    file.path(lib, "src", dll)
  ))

  for (path in candidates) {
    if (file.exists(path)) {
      dyn.load(path)
      return(TRUE)
    }
  }

  isTRUE(tryCatch({
    library.dynam(pkg, pkg, lib)
    TRUE
  }, error = function(e) FALSE))
}

.onAttach <- function (lib, pkg) {
		packageStartupMessage(
      sprintf(
        paste(
          "npRmpi %s",
          "Examples and guides at https://jeffreyracine.github.io/gallery/",
          'See also vignette("npRmpi_getting_started", package = "npRmpi")',
          sep = "\n"
        ),
        utils::packageDescription(pkg, fields = "Version")
      ),
      domain = NULL,
      appendLF = TRUE
    )
    if (isTRUE(getOption("npRmpi.conflicts.warn", TRUE)) &&
        ("package:np" %in% search()) &&
        !isTRUE(getOption("npRmpi.conflicts.warned", FALSE))) {
      packageStartupMessage("note: both 'npRmpi' and 'np' are attached; prefer explicit npRmpi:: calls to avoid masking ambiguity")
      options(npRmpi.conflicts.warned = TRUE)
    }
    if (isTRUE(getOption("npRmpi.conflicts.warn", TRUE)) &&
        ("package:Rmpi" %in% search()) &&
        !isTRUE(getOption("npRmpi.conflicts.warned.rmpi", FALSE))) {
      packageStartupMessage("note: both 'npRmpi' and 'Rmpi' are attached; prefer explicit npRmpi:: calls for estimator APIs")
      options(npRmpi.conflicts.warned.rmpi = TRUE)
    }
}

.onUnload <- function (lpath){
  tryCatch(.Call("C_np_release_static_buffers", PACKAGE = "npRmpi"), error = function(e) NULL)
  if (isTRUE(getOption("npRmpi.mpi.initialized", FALSE)))
    tryCatch(mpi.finalize(), error = function(e) NULL)
}

.onLoad <- function (lib, pkg) {
  if (!is.loaded("mpi_initialize") && !.npRmpi_try_dynload(lib = lib, pkg = pkg))
    stop("Fail to load npRmpi dynamic library.")
  if (!is.loaded("mpi_initialize"))
    stop("Probably npRmpi has been detached. Please quit R.")

  # On macOS with MPICH, repeated spawn/close cycles have been observed to
  # destabilize subsequent MPI collectives. Default to reusing spawned slaves
  # (keep them alive) unless the user explicitly disables it.
  if (is.null(getOption("npRmpi.reuse.slaves")) &&
      identical(Sys.info()[["sysname"]], "Darwin") &&
      !nzchar(Sys.getenv("NP_RMPI_NO_REUSE_SLAVES"))) {
    options(npRmpi.reuse.slaves = TRUE)
  }

	if (nzchar(Sys.getenv("NP_RMPI_SKIP_INIT"))) {
    options(npRmpi.mpi.initialized = FALSE)
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
    if (is.null(getOption("npRmpi.autodispatch")))
      options(npRmpi.autodispatch = FALSE)
    if (is.null(getOption("npRmpi.autodispatch.strict")))
      options(npRmpi.autodispatch.strict = TRUE)
    if (is.null(getOption("npRmpi.autodispatch.disable")))
      options(npRmpi.autodispatch.disable = FALSE)
    if (is.null(getOption("npRmpi.autodispatch.context")))
      options(npRmpi.autodispatch.context = FALSE)
    if (is.null(getOption("npRmpi.autodispatch.warned.nested")))
      options(npRmpi.autodispatch.warned.nested = FALSE)
    if (is.null(getOption("npRmpi.master.only")))
      options(npRmpi.master.only = FALSE)
    if (is.null(getOption("npRmpi.pool.active")))
      options(npRmpi.pool.active = FALSE)
    if (is.null(getOption("npRmpi.conflicts.warn")))
      options(npRmpi.conflicts.warn = TRUE)
    if (is.null(getOption("npRmpi.conflicts.warned")))
      options(npRmpi.conflicts.warned = FALSE)
	    if (is.null(getOption("npRmpi.conflicts.warned.rmpi")))
	      options(npRmpi.conflicts.warned.rmpi = FALSE)
    if (is.null(getOption("npRmpi.embedded.backend.version")))
      options(npRmpi.embedded.backend.version = .npRmpi_embedded_backend_version)

    setHook(packageEvent("np", "attach"),
            function(...) tryCatch(get(".npRmpi_after_np_attach",
                                       envir = asNamespace("npRmpi"),
                                       mode = "function",
                                       inherits = FALSE)(),
                                   error = function(e) invisible(e)),
            action = "append")
    setHook(packageEvent("np", "onLoad"),
            function(...) tryCatch(get(".npRmpi_after_np_load",
                                       envir = asNamespace("npRmpi"),
                                       mode = "function",
                                       inherits = FALSE)(),
                                   error = function(e) invisible(e)),
            action = "append")
    if ("np" %in% loadedNamespaces())
      .npRmpi_after_np_load()
	    return(invisible())
	  }

  if(.Call("mpidist",PACKAGE="npRmpi") == 2){
    auto.lamboot <- isTRUE(getOption("npRmpi.auto.lamboot", FALSE)) ||
      nzchar(Sys.getenv("NP_RMPI_AUTO_LAMBOOT"))
    lamnodes <- tryCatch(system("lamnodes", TRUE, ignore.stderr = TRUE),
                         error = function(e) character(0))
    if (auto.lamboot && length(lamnodes) == 0){
	    system("lamboot -H",ignore.stderr = TRUE)
    }
  }
	
  if(!.Call("mpi_initialize",PACKAGE="npRmpi"))
    stop("Cannot start MPI_Init(). Exit")
  options(npRmpi.mpi.initialized = TRUE)
  
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
  if (is.null(getOption("npRmpi.autodispatch")))
    options(npRmpi.autodispatch = FALSE)
  if (is.null(getOption("npRmpi.autodispatch.strict")))
    options(npRmpi.autodispatch.strict = TRUE)
  if (is.null(getOption("npRmpi.autodispatch.disable")))
    options(npRmpi.autodispatch.disable = FALSE)
  if (is.null(getOption("npRmpi.autodispatch.context")))
    options(npRmpi.autodispatch.context = FALSE)
  if (is.null(getOption("npRmpi.autodispatch.warned.nested")))
    options(npRmpi.autodispatch.warned.nested = FALSE)
  if (is.null(getOption("npRmpi.master.only")))
    options(npRmpi.master.only = FALSE)
  if (is.null(getOption("npRmpi.pool.active")))
    options(npRmpi.pool.active = FALSE)
  if (is.null(getOption("npRmpi.conflicts.warn")))
    options(npRmpi.conflicts.warn = TRUE)
  if (is.null(getOption("npRmpi.conflicts.warned")))
    options(npRmpi.conflicts.warned = FALSE)
  if (is.null(getOption("npRmpi.conflicts.warned.rmpi")))
    options(npRmpi.conflicts.warned.rmpi = FALSE)
  if (is.null(getOption("npRmpi.embedded.backend.version")))
    options(npRmpi.embedded.backend.version = .npRmpi_embedded_backend_version)

  setHook(packageEvent("np", "attach"),
          function(...) tryCatch(get(".npRmpi_after_np_attach",
                                     envir = asNamespace("npRmpi"),
                                     mode = "function",
                                     inherits = FALSE)(),
                                 error = function(e) invisible(e)),
          action = "append")
  setHook(packageEvent("np", "onLoad"),
          function(...) tryCatch(get(".npRmpi_after_np_load",
                                     envir = asNamespace("npRmpi"),
                                     mode = "function",
                                     inherits = FALSE)(),
                                 error = function(e) invisible(e)),
          action = "append")
  if ("np" %in% loadedNamespaces())
    .npRmpi_after_np_load()
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

npMacMseriesAccelerateMode <- function(value = getOption("np.macMseries.accelerate", "auto")) {
  if (isTRUE(value))
    return("on")
  if (identical(value, FALSE))
    return("off")
  if (is.character(value) && length(value) == 1L &&
      identical(unname(tolower(value)), "auto"))
    return("auto")
  stop("option 'np.macMseries.accelerate' must be TRUE, FALSE, or \"auto\"", call. = FALSE)
}

.npMacMseriesAccelerateOptionValue <- function(value = getOption("np.macMseries.accelerate", "auto")) {
  mode <- npMacMseriesAccelerateMode(value)
  switch(mode,
         on = TRUE,
         off = FALSE,
         auto = "auto")
}

.npRmpi_sync_runtime_options_to_slaves <- function(comm = 1L, force = FALSE) {
  if (isTRUE(getOption("npRmpi.runtime.option.sync.active", FALSE)))
    return(invisible(FALSE))

  rank <- tryCatch(as.integer(mpi.comm.rank(comm)), error = function(e) NA_integer_)
  size <- tryCatch(as.integer(mpi.comm.size(comm)), error = function(e) NA_integer_)
  if (is.na(rank) || is.na(size) || rank != 0L || size < 2L)
    return(invisible(FALSE))

  accel <- .npMacMseriesAccelerateOptionValue()
  last <- getOption("npRmpi.macMseries.accelerate.last.synced", NULL)
  if (!isTRUE(force) && identical(accel, last))
    return(invisible(FALSE))

  old <- getOption("npRmpi.runtime.option.sync.active", FALSE)
  options(npRmpi.runtime.option.sync.active = TRUE)
  on.exit(options(npRmpi.runtime.option.sync.active = old), add = TRUE)

  mpi.bcast.cmd(options,
                np.macMseries.accelerate = accel,
                rank = 0L,
                comm = comm)
  options(npRmpi.macMseries.accelerate.last.synced = accel)
  invisible(TRUE)
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
