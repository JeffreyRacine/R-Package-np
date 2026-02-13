#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: Rscript benchmarks/compare_np_serial_microbench_synth.R <pre_dir> <post_dir>")
}

pre_dir <- args[[1]]
post_dir <- args[[2]]

pre_file <- file.path(pre_dir, "synth_bench.csv")
post_file <- file.path(post_dir, "synth_bench.csv")

if (!file.exists(pre_file)) stop("Missing file: ", pre_file)
if (!file.exists(post_file)) stop("Missing file: ", post_file)

pre <- read.csv(pre_file, stringsAsFactors = FALSE)
post <- read.csv(post_file, stringsAsFactors = FALSE)

by_keys <- c("n", "phase", "task")
m <- merge(pre, post, by = by_keys, suffixes = c("_pre", "_post"))

m$mean_speedup <- m$mean_ms_pre / m$mean_ms_post
m$median_speedup <- m$median_ms_pre / m$median_ms_post
m$mean_reduction_pct <- 100.0 * (m$mean_ms_pre - m$mean_ms_post) / m$mean_ms_pre
m$median_reduction_pct <- 100.0 * (m$median_ms_pre - m$median_ms_post) / m$median_ms_pre

out <- m[, c(
  "n", "phase", "task",
  "mean_ms_pre", "mean_ms_post", "mean_speedup", "mean_reduction_pct",
  "median_ms_pre", "median_ms_post", "median_speedup", "median_reduction_pct"
)]

out <- out[order(out$n, out$phase, out$task), ]
print(out, row.names = FALSE)

out_file <- file.path(post_dir, "comparison_vs_pre_synth.csv")
write.csv(out, out_file, row.names = FALSE)
message("Wrote comparison: ", out_file)
