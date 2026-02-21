.np_warn_npplot_deprecated_once <- local({
  warned <- FALSE
  function() {
    if (warned)
      return(invisible(FALSE))
    warned <<- TRUE
    warning(
      "npplot() is a compatibility API and may be deprecated in a future release; ",
      "prefer plot(...) S3 methods.",
      call. = FALSE
    )
    invisible(TRUE)
  }
})

npplot <- function(bws = stop("'bws' has not been set"), ..., random.seed = 42,
                   .npplot.internal = FALSE){
  if (!isTRUE(.npplot.internal))
    .np_warn_npplot_deprecated_once()

  seed.state <- .np_seed_enter(random.seed)
  on.exit(.np_seed_exit(seed.state), add = TRUE)
  UseMethod("npplot", bws)
}
