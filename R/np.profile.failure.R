# Called only by the standalone profile launcher, never by shared worker loops.
.npRmpi_profile_abort_error <- function(e = NULL) {
  if (!is.null(e)) {
    text <- if (inherits(e, "condition")) conditionMessage(e) else as.character(e)
    cat(paste(text, collapse = "\n"), "\n", file = stderr())
    flush(stderr())
  }
  returned <- .Call("np_profile_job_abort", PACKAGE = "npRmpi")
  if (isTRUE(returned)) {
    # Abort was requested but returned. The job is terminal: do not attempt
    # collective finalization again or recurse through the top-level handler.
    options(npRmpi.mpi.initialized = FALSE, error = NULL)
    cat("MPI_Abort returned after a terminal profile error; exiting without MPI cleanup.\n",
        file = stderr())
    quit(save = "no", status = 1L, runLast = FALSE)
  }
  invisible(FALSE)
}

.npRmpi_profile_install_error_handler <- function() {
  owned <- .Call("np_profile_job_register", PACKAGE = "npRmpi")
  if (isTRUE(owned)) {
    # This hook runs only after ordinary R condition handlers decline an error.
    # A batch profile owns its job's terminal error policy; caught errors retain
    # their usual semantics. Session/attach and one-rank R never install it.
    options(error = function() .npRmpi_profile_abort_error(geterrmessage()))
  }
  invisible(owned)
}
