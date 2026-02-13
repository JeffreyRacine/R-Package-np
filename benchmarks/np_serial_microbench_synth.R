#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(np)
  library(microbenchmark)
})

options(np.messages = FALSE)

seed <- as.integer(Sys.getenv("NP_SYNTH_SEED", "42"))
n_values <- as.integer(strsplit(Sys.getenv("NP_SYNTH_N", "500,1000,2000"), ",")[[1]])
times_combined <- as.integer(Sys.getenv("NP_SYNTH_TIMES_COMBINED", "20"))
times_split <- as.integer(Sys.getenv("NP_SYNTH_TIMES_SPLIT", "10"))
warmup <- as.integer(Sys.getenv("NP_SYNTH_WARMUP", "1"))
mode <- Sys.getenv("NP_SYNTH_MODE", "combined") # combined|split|both
out_dir <- Sys.getenv("NP_SYNTH_OUTDIR", unset = file.path(tempdir(), "np_serial_microbench_synth"))
tasks_env <- Sys.getenv("NP_SYNTH_TASKS", "npreg,npcdens,npcdist,npudens,npudist")
tasks <- trimws(strsplit(tasks_env, ",")[[1]])
valid_tasks <- c("npreg", "npcdens", "npcdist", "npudens", "npudist")

if (length(n_values) < 1L || any(is.na(n_values)) || any(n_values <= 10L)) {
  stop("NP_SYNTH_N must contain positive sample sizes > 10")
}
if (times_combined < 1L || times_split < 1L || warmup < 0L) {
  stop("Invalid settings: times must be >= 1 and warmup >= 0")
}
if (!(mode %in% c("combined", "split", "both"))) {
  stop("NP_SYNTH_MODE must be one of: combined, split, both")
}
if (length(tasks) < 1L || any(!(tasks %in% valid_tasks))) {
  stop("NP_SYNTH_TASKS must be a comma-separated subset of: ",
       paste(valid_tasks, collapse = ","))
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mk_data <- function(n, seed_offset = 0L) {
  set.seed(seed + seed_offset + n)
  x1 <- runif(n)
  z1 <- rbinom(n, 1L, .5)
  z2 <- rbinom(n, 1L, .5)
  y <- cos(2 * pi * x1) + z1 + rnorm(n, sd = .25)
  data.frame(
    y = y,
    x1 = x1,
    z1 = factor(z1),
    z1o = ordered(z1),
    z2 = factor(z2)
  )
}

summarize_bm <- function(bm, n, phase, task) {
  sm <- summary(bm)
  data.frame(
    n = n,
    phase = phase,
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

bench_combined <- function(dat, n) {
  exprs_all <- list(
    npreg_combined = function() npreg(y ~ x1 + z1 + z2, data = dat, nmulti = 1),
    npcdens_combined = function() npcdens(y ~ x1 + z1 + z2, data = dat, nmulti = 1),
    npcdist_combined = function() npcdist(y ~ x1 + z1 + z2, data = dat, nmulti = 1),
    npudens_combined = function() npudens(~ y + x1 + z1o, data = dat, nmulti = 1),
    npudist_combined = function() npudist(~ y + x1 + z1o, data = dat, nmulti = 1)
  )
  exprs <- exprs_all[paste0(tasks, "_combined")]

  out <- list()
  for (nm in names(exprs)) {
    run_warmup(exprs[[nm]], warmup)
    bm <- microbenchmark(exprs[[nm]](), times = times_combined, unit = "ms")
    out[[nm]] <- summarize_bm(bm, n, "combined", nm)
  }
  do.call(rbind, out)
}

bench_split <- function(dat, n) {
  bw_all <- list(
    npreg_bw = function() npregbw(y ~ x1 + z1 + z2, data = dat, nmulti = 1),
    npcdens_bw = function() npcdensbw(y ~ x1 + z1 + z2, data = dat, nmulti = 1),
    npcdist_bw = function() npcdistbw(y ~ x1 + z1 + z2, data = dat, nmulti = 1),
    npudens_bw = function() npudensbw(~ y + x1 + z1o, data = dat, nmulti = 1),
    npudist_bw = function() npudistbw(~ y + x1 + z1o, data = dat, nmulti = 1)
  )
  bw_exprs <- bw_all[paste0(tasks, "_bw")]

  fit_all <- list(
    npreg_fit = function() {
      bw <- npregbw(y ~ x1 + z1 + z2, data = dat, nmulti = 1)
      function() npreg(bws = bw)
    },
    npcdens_fit = function() {
      bw <- npcdensbw(y ~ x1 + z1 + z2, data = dat, nmulti = 1)
      function() npcdens(bws = bw)
    },
    npcdist_fit = function() {
      bw <- npcdistbw(y ~ x1 + z1 + z2, data = dat, nmulti = 1)
      function() npcdist(bws = bw)
    },
    npudens_fit = function() {
      bw <- npudensbw(~ y + x1 + z1o, data = dat, nmulti = 1)
      function() npudens(bws = bw)
    },
    npudist_fit = function() {
      bw <- npudistbw(~ y + x1 + z1o, data = dat, nmulti = 1)
      function() npudist(bws = bw)
    }
  )
  fit_prep <- fit_all[paste0(tasks, "_fit")]

  out <- list()

  for (nm in names(bw_exprs)) {
    run_warmup(bw_exprs[[nm]], warmup)
    bm <- microbenchmark(bw_exprs[[nm]](), times = times_split, unit = "ms")
    out[[nm]] <- summarize_bm(bm, n, "split_bw", nm)
  }

  for (nm in names(fit_prep)) {
    fit_fun <- fit_prep[[nm]]()
    run_warmup(fit_fun, warmup)
    bm <- microbenchmark(fit_fun(), times = times_split, unit = "ms")
    out[[nm]] <- summarize_bm(bm, n, "split_fit", nm)
  }

  do.call(rbind, out)
}

all_results <- list()
for (n in n_values) {
  message("Running synthetic benchmarks for n = ", n)
  dat <- mk_data(n)
  if (mode %in% c("combined", "both")) {
    all_results[[paste0("combined_", n)]] <- bench_combined(dat, n)
  }
  if (mode %in% c("split", "both")) {
    all_results[[paste0("split_", n)]] <- bench_split(dat, n)
  }
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
