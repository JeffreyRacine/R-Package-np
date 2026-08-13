locate_rbandwidth_plot_engine <- function() {
  candidates <- c(
    file.path("R", "np.plot.engine.rbandwidth.R"),
    file.path("..", "..", "R", "np.plot.engine.rbandwidth.R")
  )
  hit <- candidates[file.exists(candidates)]
  if (!length(hit))
    return(NA_character_)
  normalizePath(hit[[1L]], mustWork = TRUE)
}

test_that("every direct rbandwidth asymptotic fit requests standard errors", {
  path <- locate_rbandwidth_plot_engine()
  skip_if(is.na(path), "source tree unavailable")
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")

  direct.calls <- gregexpr("npreg(", text, fixed = TRUE)[[1L]]
  direct.count <- if (identical(direct.calls, -1L)) 0L else length(direct.calls)
  se.requests <- gregexpr("se = TRUE,", text, fixed = TRUE)[[1L]]
  se.count <- if (identical(se.requests, -1L)) 0L else length(se.requests)

  expect_identical(direct.count, 3L)
  expect_identical(se.count, direct.count)
})
