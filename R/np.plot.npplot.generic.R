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
  .np_plot_call_method(.np_plot_compat_dispatch, bws = bws, ...,
                       random.seed = random.seed, where = "npplot()")
}
