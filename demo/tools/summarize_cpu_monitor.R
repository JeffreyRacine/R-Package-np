#!/usr/bin/env Rscript

cmd <- commandArgs(FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)
if (!length(file_arg)) stop("cannot locate wrapper path", call. = FALSE)
wrapper <- normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(wrapper), "..", ".."), mustWork = TRUE)
source(file.path(repo_root, "inst", "demo_tools", "summarize_cpu_monitor.R"))
