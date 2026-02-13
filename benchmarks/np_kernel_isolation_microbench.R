#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(np)
  library(microbenchmark)
})

options(np.messages = FALSE)

seed <- as.integer(Sys.getenv("NP_ISO_SEED", "20260212"))
n_values <- as.integer(strsplit(Sys.getenv("NP_ISO_N", "500,1000,2000"), ",")[[1]])
times_combined <- as.integer(Sys.getenv("NP_ISO_TIMES", "10"))
warmup <- as.integer(Sys.getenv("NP_ISO_WARMUP", "2"))
out_dir <- Sys.getenv("NP_ISO_OUTDIR", unset = file.path(tempdir(), "np_kernel_isolation"))
include_x <- tolower(Sys.getenv("NP_ISO_INCLUDE_X", "false")) %in% c("1", "true", "yes")

if (length(n_values) < 1L || any(is.na(n_values)) || any(n_values <= 10L)) {
  stop("NP_ISO_N must contain positive sample sizes > 10")
}
if (times_combined < 1L || warmup < 0L) {
  stop("Invalid settings: times must be >= 1 and warmup >= 0")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mk_data <- function(n, seed_offset = 0L, include_x = FALSE) {
  set.seed(seed + seed_offset + n)

  z1 <- rbinom(n, 3L, .5)
  z2 <- rbinom(n, 3L, .5)

  if (include_x) {
    x1 <- runif(n)
    y <- cos(2 * pi * x1) + 0.5 * z1 + 0.25 * z2 + rnorm(n, sd = .25)
    data.frame(
      y = y,
      x1 = x1,
      zf1 = factor(z1),
      zf2 = factor(z2),
      zo1 = ordered(z1),
      zo2 = ordered(z2)
    )
  } else {
    y <- sin(z1) + 0.3 * z2 + rnorm(n, sd = .25)
    data.frame(
      y = y,
      zf1 = factor(z1),
      zf2 = factor(z2),
      zo1 = ordered(z1),
      zo2 = ordered(z2)
    )
  }
}

summarize_bm <- function(bm, n, task) {
  sm <- summary(bm)
  data.frame(
    n = n,
    phase = "combined",
    task = task,
    min_ms = sm$min,
    lq_ms = sm$lq,
    mean_ms = sm$mean,
    median_ms = sm$median,
    uq_ms = sm$uq,
    max_ms = sm$max,
    stringsAsFactors = FALSE
  )
}

run_warmup <- function(expr_fun, warmup_n) {
  if (warmup_n < 1L) return(invisible(NULL))
  for (i in seq_len(warmup_n)) expr_fun()
  invisible(NULL)
}

bench_one_n <- function(dat, n, include_x = FALSE) {
  if (include_x) {
    exprs <- list(
      npreg_combined = function() npreg(y ~ x1 + zf1 + zo1, data = dat, nmulti = 1),
      npcdens_combined = function() npcdens(y ~ x1 + zf1 + zo1, data = dat, nmulti = 1),
      npudens_combined = function() npudens(~ y + x1 + zo1, data = dat, nmulti = 1),
      npudist_combined = function() npudist(~ y + x1 + zo1, data = dat, nmulti = 1)
    )
  } else {
    exprs <- list(
      npreg_combined = function() npreg(y ~ zf1 + zo1 + zf2 + zo2, data = dat, nmulti = 1),
      npcdens_combined = function() npcdens(y ~ zf1 + zo1 + zf2 + zo2, data = dat, nmulti = 1),
      npudens_combined = function() npudens(~ y + zo1 + zo2, data = dat, nmulti = 1),
      npudist_combined = function() npudist(~ y + zo1 + zo2, data = dat, nmulti = 1)
    )
  }

  out <- list()
  for (nm in names(exprs)) {
    run_warmup(exprs[[nm]], warmup)
    bm <- microbenchmark(exprs[[nm]](), times = times_combined, unit = "ms")
    out[[nm]] <- summarize_bm(bm, n, nm)
  }
  do.call(rbind, out)
}

all_results <- list()
for (n in n_values) {
  message("Running kernel-isolation synthetic benchmarks for n = ", n)
  dat <- mk_data(n, include_x = include_x)
  all_results[[as.character(n)]] <- bench_one_n(dat, n, include_x = include_x)
}

res <- do.call(rbind, all_results)
res <- res[order(res$n, res$phase, res$task), ]

csv_file <- file.path(out_dir, "synth_bench.csv")
meta_file <- file.path(out_dir, "session_info.txt")

write.csv(res, csv_file, row.names = FALSE)
writeLines(capture.output(sessionInfo()), con = meta_file)

print(res[, c("n", "phase", "task", "mean_ms", "median_ms", "min_ms", "max_ms")], row.names = FALSE)
message("")
message("Outputs:")
message("  ", csv_file)
message("  ", meta_file)
