.npRmpi_embedded_backend_version <- "0.7-3.3"

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
      paste0(
        "Parallel Nonparametric Kernel Methods for Mixed Datatypes (version 0.70-1)",
        " + embedded Rmpi backend ",
        .npRmpi_embedded_backend_version,
        "\n[vignette(\"np_faq\",package=\"npRmpi\") provides answers to frequently asked questions]",
        "\n[vignette(\"npRmpi\",package=\"npRmpi\") an overview]",
        "\n[vignette(\"entropy_np\",package=\"npRmpi\") an overview of entropy-based methods]"
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
  library.dynam.unload("npRmpi", libpath=lpath) 
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
      options(np.tree = FALSE)
    if (is.null(getOption("np.largeh.rel.tol")))
      options(np.largeh.rel.tol = 1e-3)
    if (is.null(getOption("np.disc.upper.rel.tol")))
      options(np.disc.upper.rel.tol = 1e-2)
    if (is.null(getOption("np.groupcv.fast")))
      options(np.groupcv.fast = TRUE)
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
    if (is.null(getOption("npRmpi.conflicts.warn")))
      options(npRmpi.conflicts.warn = TRUE)
    if (is.null(getOption("npRmpi.conflicts.warned")))
      options(npRmpi.conflicts.warned = FALSE)
	    if (is.null(getOption("npRmpi.conflicts.warned.rmpi")))
	      options(npRmpi.conflicts.warned.rmpi = FALSE)
    if (is.null(getOption("npRmpi.embedded.backend.version")))
      options(npRmpi.embedded.backend.version = .npRmpi_embedded_backend_version)

    setHook(packageEvent("np", "attach"),
            function(...) tryCatch(get(".npRmpi_warn_pkg_conflict_once",
                                       envir = asNamespace("npRmpi"),
                                       mode = "function",
                                       inherits = FALSE)(),
                                   error = function(e) invisible(e)),
            action = "append")
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
    options(np.tree = FALSE)
  if (is.null(getOption("np.largeh.rel.tol")))
    options(np.largeh.rel.tol = 1e-3)
  if (is.null(getOption("np.disc.upper.rel.tol")))
    options(np.disc.upper.rel.tol = 1e-2)
  if (is.null(getOption("np.groupcv.fast")))
    options(np.groupcv.fast = TRUE)
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
  if (is.null(getOption("npRmpi.conflicts.warn")))
    options(npRmpi.conflicts.warn = TRUE)
  if (is.null(getOption("npRmpi.conflicts.warned")))
    options(npRmpi.conflicts.warned = FALSE)
  if (is.null(getOption("npRmpi.conflicts.warned.rmpi")))
    options(npRmpi.conflicts.warned.rmpi = FALSE)
  if (is.null(getOption("npRmpi.embedded.backend.version")))
    options(npRmpi.embedded.backend.version = .npRmpi_embedded_backend_version)

  setHook(packageEvent("np", "attach"),
          function(...) tryCatch(get(".npRmpi_warn_pkg_conflict_once",
                                     envir = asNamespace("npRmpi"),
                                     mode = "function",
                                     inherits = FALSE)(),
                                 error = function(e) invisible(e)),
          action = "append")
}
