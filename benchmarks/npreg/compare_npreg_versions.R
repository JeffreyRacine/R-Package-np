#!/usr/bin/env Rscript

parse_args <- function(args) {
  out <- list(
    raw_a = "",
    raw_b = "",
    label_a = "A",
    label_b = "B",
    out_timing = "/tmp/npreg_version_timing_compare.csv",
    out_objective = "/tmp/npreg_version_objective_compare.csv"
  )
  for (a in args) {
    if (!grepl("^--", a)) stop("Bad arg: ", a)
    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
    key <- kv[1]
    val <- if (length(kv) > 1L) kv[2] else ""
    if (key %in% names(out)) out[[key]] <- val else stop("Unknown arg: ", key)
  }
  if (!nzchar(out$raw_a) || !nzchar(out$raw_b)) stop("--raw_a and --raw_b are required")
  out
}

mean_or_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
median_or_na <- function(x) if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)

mk_timing <- function(df, label) {
  data.frame(
    function_name = c("npregbw", "npreg", "npreg_total"),
    mean = c(mean_or_na(df$elapsed_bw), mean_or_na(df$elapsed_fit), mean_or_na(df$elapsed_total)),
    median = c(median_or_na(df$elapsed_bw), median_or_na(df$elapsed_fit), median_or_na(df$elapsed_total)),
    label = label,
    stringsAsFactors = FALSE
  )
}

pct_faster <- function(a, b) {
  if (is.na(a) || is.na(b) || a == 0) return(NA_real_)
  100 * (a - b) / a
}

main <- function() {
  cfg <- parse_args(commandArgs(trailingOnly = TRUE))
  a <- read.csv(cfg$raw_a, stringsAsFactors = FALSE)
  b <- read.csv(cfg$raw_b, stringsAsFactors = FALSE)

  ta <- mk_timing(a, cfg$label_a)
  tb <- mk_timing(b, cfg$label_b)
  tm <- merge(ta, tb, by = "function_name", suffixes = c("_a", "_b"))
  timing <- data.frame(
    function_name = tm$function_name,
    label_a = cfg$label_a,
    mean_a = tm$mean_a,
    median_a = tm$median_a,
    label_b = cfg$label_b,
    mean_b = tm$mean_b,
    median_b = tm$median_b,
    pct_faster_mean_b_vs_a = mapply(pct_faster, tm$mean_a, tm$mean_b),
    pct_faster_median_b_vs_a = mapply(pct_faster, tm$median_a, tm$median_b),
    stringsAsFactors = FALSE
  )
  write.csv(timing, cfg$out_timing, row.names = FALSE)

  keys <- intersect(names(a), names(b))
  keys <- keys[keys %in% c("n","regtype","bwmethod","nmulti","ckertype","np_tree","seed_policy","seed","iter")]
  m <- merge(a, b, by = keys, suffixes = c("_a", "_b"), all = TRUE)

  absdiff <- function(x, y) ifelse(is.na(x) | is.na(y), NA_real_, abs(x - y))
  m$diff_fval <- absdiff(m$fval_a, m$fval_b)
  m$diff_ifval <- absdiff(m$ifval_a, m$ifval_b)
  m$diff_num_fval <- absdiff(m$num_fval_a, m$num_fval_b)
  m$same_bw <- m$bw_a == m$bw_b
  m$same_ok <- m$ok_a == m$ok_b

  objective <- data.frame(
    label_a = cfg$label_a,
    label_b = cfg$label_b,
    rows_compared = nrow(m),
    ok_match_rate = mean(m$same_ok, na.rm = TRUE),
    bw_match_rate = mean(m$same_bw, na.rm = TRUE),
    max_abs_diff_fval = if (all(is.na(m$diff_fval))) NA_real_ else max(m$diff_fval, na.rm = TRUE),
    max_abs_diff_ifval = if (all(is.na(m$diff_ifval))) NA_real_ else max(m$diff_ifval, na.rm = TRUE),
    max_abs_diff_num_fval = if (all(is.na(m$diff_num_fval))) NA_real_ else max(m$diff_num_fval, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  write.csv(objective, cfg$out_objective, row.names = FALSE)

  cat("timing_csv=", cfg$out_timing, "\n", sep = "")
  cat("objective_csv=", cfg$out_objective, "\n", sep = "")
}

main()
