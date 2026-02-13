#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(np)
  library(microbenchmark)
})

options(np.messages = FALSE)

seed <- as.integer(Sys.getenv("NP_BENCH_SEED", "42"))
times_fit <- as.integer(Sys.getenv("NP_BENCH_TIMES_FIT", "100"))
times_bw <- as.integer(Sys.getenv("NP_BENCH_TIMES_BW", "20"))
warmup <- as.integer(Sys.getenv("NP_BENCH_WARMUP", "2"))
tol <- as.numeric(Sys.getenv("NP_BENCH_TOL", "1e-10"))
out_dir <- Sys.getenv("NP_BENCH_OUTDIR", unset = file.path(tempdir(), "np_serial_microbench"))

if (times_fit < 1L || times_bw < 1L || warmup < 0L) {
  stop("Invalid benchmark settings: require times_fit >= 1, times_bw >= 1, warmup >= 0")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)

# Data from package examples/man pages.
data("wage1", package = "np")
data("Italy", package = "np")
data("oecdpanel", package = "np")

reg_dat <- within(wage1[seq_len(min(220L, nrow(wage1))), c("lwage", "educ", "exper", "tenure", "female", "married")], {
  female <- factor(female)
  married <- factor(married)
})

dens_dat <- Italy[seq_len(min(380L, nrow(Italy))), c("year", "gdp")]
dens_dat$year <- ordered(dens_dat$year)

cond_dat <- within(oecdpanel[seq_len(min(260L, nrow(oecdpanel))), c("growth", "oecd", "initgdp", "popgro")], {
  oecd <- factor(oecd)
})

message("Preparing bandwidth objects (single-start for repeatability/speed)...")
bw_npreg <- npregbw(lwage ~ educ + exper + tenure + female + married, data = reg_dat, nmulti = 1)
bw_udens <- npudensbw(formula = ~ordered(year) + gdp, data = dens_dat, nmulti = 1)
bw_udist <- npudistbw(formula = ~ordered(year) + gdp, data = dens_dat, nmulti = 1)
bw_cdens <- npcdensbw(formula = growth ~ oecd + initgdp + popgro, data = cond_dat, nmulti = 1)
bw_cdist <- npcdistbw(formula = growth ~ oecd + initgdp + popgro, data = cond_dat, nmulti = 1)

fit_npreg_ref <- fitted(npreg(bws = bw_npreg, gradients = FALSE))
fit_udens_ref <- fitted(npudens(bws = bw_udens))
fit_udist_ref <- fitted(npudist(bws = bw_udist))
fit_cdens_ref <- fitted(npcdens(bws = bw_cdens))
fit_cdist_ref <- fitted(npcdist(bws = bw_cdist))

validate <- function(name, x, y, tolerance = tol) {
  ae <- all.equal(x, y, tolerance = tolerance, check.attributes = FALSE)
  ok <- isTRUE(ae)
  data.frame(task = name, ok = ok, detail = if (ok) "OK" else as.character(ae), stringsAsFactors = FALSE)
}

run_warmup <- function(fun, n) {
  if (n < 1L) return(invisible(NULL))
  for (i in seq_len(n)) fun()
  invisible(NULL)
}

bench_and_summarize <- function(name, expr, times) {
  bm <- microbenchmark(eval.parent(substitute(expr)), times = times, unit = "ms")
  sm <- summary(bm)
  data.frame(
    task = name,
    n = nrow(sm),
    min_ms = sm$min,
    lq_ms = sm$lq,
    mean_ms = sm$mean,
    median_ms = sm$median,
    uq_ms = sm$uq,
    max_ms = sm$max,
    stringsAsFactors = FALSE
  )
}

message("Warm-up...")
run_warmup(function() npreg(bws = bw_npreg, gradients = FALSE), warmup)
run_warmup(function() npudens(bws = bw_udens), warmup)
run_warmup(function() npudist(bws = bw_udist), warmup)
run_warmup(function() npcdens(bws = bw_cdens), warmup)
run_warmup(function() npcdist(bws = bw_cdist), warmup)

message("Running fit-path benchmarks...")
fit_results <- do.call(
  rbind,
  list(
    bench_and_summarize("npreg_fit", npreg(bws = bw_npreg, gradients = FALSE), times_fit),
    bench_and_summarize("npudens_fit", npudens(bws = bw_udens), times_fit),
    bench_and_summarize("npudist_fit", npudist(bws = bw_udist), times_fit),
    bench_and_summarize("npcdens_fit", npcdens(bws = bw_cdens), times_fit),
    bench_and_summarize("npcdist_fit", npcdist(bws = bw_cdist), times_fit)
  )
)

message("Running bandwidth-search benchmarks...")
bw_results <- do.call(
  rbind,
  list(
    bench_and_summarize("npreg_bw", npregbw(lwage ~ educ + exper + tenure + female + married, data = reg_dat, nmulti = 1), times_bw),
    bench_and_summarize("npudens_bw", npudensbw(formula = ~ordered(year) + gdp, data = dens_dat, nmulti = 1), times_bw),
    bench_and_summarize("npudist_bw", npudistbw(formula = ~ordered(year) + gdp, data = dens_dat, nmulti = 1), times_bw),
    bench_and_summarize("npcdens_bw", npcdensbw(formula = growth ~ oecd + initgdp + popgro, data = cond_dat, nmulti = 1), times_bw),
    bench_and_summarize("npcdist_bw", npcdistbw(formula = growth ~ oecd + initgdp + popgro, data = cond_dat, nmulti = 1), times_bw)
  )
)

message("Validating deterministic equivalence under fixed bandwidths...")
validation <- do.call(
  rbind,
  list(
    validate("npreg_fit", fit_npreg_ref, fitted(npreg(bws = bw_npreg, gradients = FALSE))),
    validate("npudens_fit", fit_udens_ref, fitted(npudens(bws = bw_udens))),
    validate("npudist_fit", fit_udist_ref, fitted(npudist(bws = bw_udist))),
    validate("npcdens_fit", fit_cdens_ref, fitted(npcdens(bws = bw_cdens))),
    validate("npcdist_fit", fit_cdist_ref, fitted(npcdist(bws = bw_cdist)))
  )
)

fit_file <- file.path(out_dir, "fit_bench.csv")
bw_file <- file.path(out_dir, "bw_bench.csv")
validation_file <- file.path(out_dir, "validation.csv")
meta_file <- file.path(out_dir, "session_info.txt")

write.csv(fit_results, fit_file, row.names = FALSE)
write.csv(bw_results, bw_file, row.names = FALSE)
write.csv(validation, validation_file, row.names = FALSE)
writeLines(capture.output(sessionInfo()), con = meta_file)

message("")
message("Fit-path summary (ms):")
print(fit_results[, c("task", "mean_ms", "median_ms", "min_ms", "max_ms")], row.names = FALSE)
message("")
message("Bandwidth-search summary (ms):")
print(bw_results[, c("task", "mean_ms", "median_ms", "min_ms", "max_ms")], row.names = FALSE)
message("")
message("Validation:")
print(validation, row.names = FALSE)
message("")
message("Outputs:")
message("  ", fit_file)
message("  ", bw_file)
message("  ", validation_file)
message("  ", meta_file)
