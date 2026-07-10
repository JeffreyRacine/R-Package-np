#!/usr/bin/env Rscript

current_script_file <- function() {
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    ofile <- frames[[i]]$ofile
    if (!is.null(ofile) && length(ofile) == 1L && nzchar(ofile)) {
      return(normalizePath(ofile, mustWork = TRUE))
    }
  }

  file_arg <- "--file="
  args <- commandArgs(trailingOnly = FALSE)
  hit <- args[startsWith(args, file_arg)]
  if (length(hit) > 0L) {
    path <- sub(file_arg, "", hit[[1L]], fixed = TRUE)
    if (nzchar(path) && !identical(path, "-") && file.exists(path)) {
      return(normalizePath(path, mustWork = TRUE))
    }
  }

  stop("Cannot determine legacy one-off benchmark script path")
}

canonical <- file.path(
  dirname(dirname(current_script_file())),
  "perf",
  "oneoff",
  "bench_oneoff_param_nprmpi.R"
)

source(canonical, local = TRUE)

if (sys.nframe() == 0) main()
