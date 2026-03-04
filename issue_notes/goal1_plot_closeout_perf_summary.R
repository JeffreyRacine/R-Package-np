#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  i <- match(name, args)
  if (is.na(i) || i == length(args)) return(default)
  args[[i + 1L]]
}

baseline_csv <- get_arg("--baseline-csv", "")
candidate_csv <- get_arg("--candidate-csv", "")
out_dir <- get_arg("--out-dir", tempfile("goal1_perf_summary_"))

if (!nzchar(baseline_csv) || !nzchar(candidate_csv))
  stop("both --baseline-csv and --candidate-csv are required")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

base <- read.csv(baseline_csv, stringsAsFactors = FALSE)
cand <- read.csv(candidate_csv, stringsAsFactors = FALSE)

bad_base <- subset(base, status != "ok" | nonfinite)
bad_cand <- subset(cand, status != "ok" | nonfinite)

if (nrow(bad_base) > 0L || nrow(bad_cand) > 0L) {
  verdict <- "FAIL"
  reason <- sprintf("bad rows baseline=%d candidate=%d", nrow(bad_base), nrow(bad_cand))
  writeLines(
    c(
      sprintf("verdict=%s", verdict),
      sprintf("reason=%s", reason),
      sprintf("baseline_csv=%s", baseline_csv),
      sprintf("candidate_csv=%s", candidate_csv)
    ),
    con = file.path(out_dir, "perf_verdict.txt")
  )
  quit(save = "no", status = 1L)
}

keep_cols <- c("scenario", "seed", "elapsed_sec")
base2 <- base[, keep_cols]
cand2 <- cand[, keep_cols]
names(base2)[3L] <- "elapsed_base"
names(cand2)[3L] <- "elapsed_cand"

merged <- merge(base2, cand2, by = c("scenario", "seed"), all = FALSE)
merged$delta_sec <- merged$elapsed_cand - merged$elapsed_base
merged$delta_rel <- merged$delta_sec / pmax(merged$elapsed_base, .Machine$double.eps)

boot_median_ci <- function(x, B = 2000L, conf = 0.95) {
  if (length(x) < 2L) return(c(NA_real_, NA_real_))
  set.seed(271828L)
  med <- replicate(B, median(sample(x, replace = TRUE)))
  alpha <- (1 - conf) / 2
  as.numeric(quantile(med, probs = c(alpha, 1 - alpha), names = FALSE))
}

summarize_one <- function(df) {
  d <- df$delta_sec
  eb <- df$elapsed_base
  ec <- df$elapsed_cand
  mean_ci <- if (length(d) >= 2L) t.test(d)$conf.int else c(NA_real_, NA_real_)
  med_ci <- boot_median_ci(d)

  base_mean <- mean(eb)
  cand_mean <- mean(ec)
  base_median <- median(eb)
  cand_median <- median(ec)
  base_p95 <- as.numeric(quantile(eb, 0.95, names = FALSE))
  cand_p95 <- as.numeric(quantile(ec, 0.95, names = FALSE))

  mean_delta <- mean(d)
  median_delta <- median(d)
  p90_delta <- as.numeric(quantile(d, 0.90, names = FALSE))
  p95_delta <- as.numeric(quantile(d, 0.95, names = FALSE))
  sd_delta <- sd(d)
  iqr_delta <- IQR(d)
  mad_delta <- mad(d, center = median(d))

  mean_rel <- mean_delta / pmax(base_mean, .Machine$double.eps)
  median_rel <- median_delta / pmax(base_median, .Machine$double.eps)
  p95_rel <- (cand_p95 - base_p95) / pmax(base_p95, .Machine$double.eps)

  center_abs_mei <- 0.10
  center_rel_mei <- 0.10
  tail_abs_mei <- 0.20
  tail_rel_mei <- 0.15

  regressed <- (mean_delta > center_abs_mei && mean_rel > center_rel_mei) ||
    (median_delta > center_abs_mei && median_rel > center_rel_mei) ||
    ((cand_p95 - base_p95) > tail_abs_mei && p95_rel > tail_rel_mei)

  data.frame(
    scenario = df$scenario[1L],
    n_pairs = nrow(df),
    base_mean = base_mean,
    cand_mean = cand_mean,
    base_median = base_median,
    cand_median = cand_median,
    mean_delta = mean_delta,
    median_delta = median_delta,
    p90_delta = p90_delta,
    p95_delta = p95_delta,
    sd_delta = sd_delta,
    iqr_delta = iqr_delta,
    mad_delta = mad_delta,
    mean_ci_lo = as.numeric(mean_ci[1L]),
    mean_ci_hi = as.numeric(mean_ci[2L]),
    median_ci_lo = med_ci[1L],
    median_ci_hi = med_ci[2L],
    mean_rel = mean_rel,
    median_rel = median_rel,
    p95_rel = p95_rel,
    regressed = regressed,
    stringsAsFactors = FALSE
  )
}

spl <- split(merged, merged$scenario)
sum_rows <- lapply(spl, summarize_one)
summary_df <- do.call(rbind, sum_rows)
summary_df <- summary_df[order(summary_df$scenario), , drop = FALSE]

summary_csv <- file.path(out_dir, "perf_summary.csv")
runs_csv <- file.path(out_dir, "perf_pairs.csv")
write.csv(summary_df, summary_csv, row.names = FALSE)
write.csv(merged, runs_csv, row.names = FALSE)

any_reg <- any(summary_df$regressed)
verdict <- if (any_reg) "FAIL" else "PASS"

writeLines(
  c(
    sprintf("verdict=%s", verdict),
    sprintf("baseline_csv=%s", baseline_csv),
    sprintf("candidate_csv=%s", candidate_csv),
    sprintf("pairs_csv=%s", runs_csv),
    sprintf("summary_csv=%s", summary_csv),
    sprintf("scenarios=%d", nrow(summary_df)),
    sprintf("regressed_scenarios=%d", sum(summary_df$regressed))
  ),
  con = file.path(out_dir, "perf_verdict.txt")
)

cat(sprintf("perf_summary_csv=%s\n", summary_csv))
cat(sprintf("perf_pairs_csv=%s\n", runs_csv))
cat(sprintf("perf_verdict=%s\n", verdict))

if (any_reg) quit(save = "no", status = 1L)
quit(save = "no", status = 0L)
