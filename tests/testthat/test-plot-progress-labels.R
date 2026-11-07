test_that("inference aborts retire the line without replaying stale preparation", {
  old <- options(np.messages=TRUE,np.plot.progress=TRUE)
  on.exit(options(old),add=TRUE)
  for(label in c("Plot inference","Plot bootstrap")) {
    z <- capture_progress_shadow_trace({
      expect_error(npRmpi:::.np_plot_progress_run({
        npRmpi:::.np_plot_progress_notify(stage="preparing",force=TRUE)
        stop("original diagnostic")
      },enabled=TRUE,label=label),"original diagnostic")
    },force_renderer="single_line")
    abort <- Filter(function(x)x$event=="abort",z$trace)
    expect_length(abort,1L)
    expect_match(abort[[1]]$line,paste(label,"interrupted"),fixed=TRUE)
    expect_false(grepl("targets pending",abort[[1]]$line,fixed=TRUE))
    expect_false(any(vapply(z$trace,function(x)x$event=="finish",TRUE)))
    expect_null(npRmpi:::.np_plot_progress_runtime$context)
    expect_null(npRmpi:::.np_progress_registry$active_id)
  }
})
