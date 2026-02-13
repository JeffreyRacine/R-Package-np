#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: Rscript benchmarks/compare_np_serial_microbench.R <pre_dir> <post_dir>")
}

pre_dir <- args[[1]]
post_dir <- args[[2]]

read_bench <- function(dir_path, file_name) {
  f <- file.path(dir_path, file_name)
  if (!file.exists(f)) stop("Missing file: ", f)
  read.csv(f, stringsAsFactors = FALSE)
}

merge_speed <- function(pre_df, post_df, label) {
  m <- merge(pre_df, post_df, by = "task", suffixes = c("_pre", "_post"))
  m$type <- label
  m$mean_speedup <- m$mean_ms_pre / m$mean_ms_post
  m$median_speedup <- m$median_ms_pre / m$median_ms_post
  m$mean_reduction_pct <- 100.0 * (m$mean_ms_pre - m$mean_ms_post) / m$mean_ms_pre
  m$median_reduction_pct <- 100.0 * (m$median_ms_pre - m$median_ms_post) / m$median_ms_pre
  m[, c(
    "type", "task",
    "mean_ms_pre", "mean_ms_post", "mean_speedup", "mean_reduction_pct",
    "median_ms_pre", "median_ms_post", "median_speedup", "median_reduction_pct"
  )]
}

fit_pre <- read_bench(pre_dir, "fit_bench.csv")
fit_post <- read_bench(post_dir, "fit_bench.csv")
bw_pre <- read_bench(pre_dir, "bw_bench.csv")
bw_post <- read_bench(post_dir, "bw_bench.csv")

fit_cmp <- merge_speed(fit_pre, fit_post, "fit")
bw_cmp <- merge_speed(bw_pre, bw_post, "bw")
cmp <- rbind(fit_cmp, bw_cmp)

print(cmp, row.names = FALSE)

out_file <- file.path(post_dir, "comparison_vs_pre.csv")
write.csv(cmp, out_file, row.names = FALSE)
message("Wrote comparison: ", out_file)
