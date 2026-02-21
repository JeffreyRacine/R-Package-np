#!/usr/bin/env Rscript

parse_args <- function(args) {
  out <- list(
    raw_a = "",
    raw_b = "",
    label_a = "A",
    label_b = "B",
    out_timing = "/tmp/npplreg_version_timing_compare.csv",
    out_objective = "/tmp/npplreg_version_objective_compare.csv",
    out_combo_timing = "/tmp/npplreg_version_combo_timing_compare.csv",
    out_combo_objective = "/tmp/npplreg_version_combo_objective_compare.csv"
  )
  for (a in args) {
    if (!grepl("^--", a)) stop("Bad arg: ", a)
    raw <- sub("^--", "", a)
    key <- sub("=.*$", "", raw)
    val <- if (grepl("=", raw, fixed = TRUE)) sub("^[^=]*=", "", raw) else ""
    if (key %in% names(out)) out[[key]] <- val else stop("Unknown arg: ", key)
  }
  if (!nzchar(out$raw_a) || !nzchar(out$raw_b)) stop("--raw_a and --raw_b are required")
  out
}

mean_or_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
median_or_na <- function(x) if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
max_or_na <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

normalize_bool <- function(x) {
  if (is.logical(x)) return(x)
  tolower(as.character(x)) %in% c("true", "t", "1")
}

prep_raw <- function(df) {
  if (!("num_fval" %in% names(df))) df$num_fval <- NA_real_
  if (!("ok" %in% names(df))) df$ok <- NA
  if (!("bw" %in% names(df))) df$bw <- ""
  if (!("np_tree" %in% names(df))) df$np_tree <- NA
  if (!("log_likelihood" %in% names(df))) df$log_likelihood <- NA_real_
  df$np_tree <- normalize_bool(df$np_tree)
  df$ok <- normalize_bool(df$ok)
  df
}

combo_keys <- function(df) {
  keep <- c("n", "bwmethod", "nmulti", "ckertype", "np_tree", "seed_policy")
  keep[keep %in% names(df)]
}

match_keys <- function(a, b) {
  keep <- c("n", "bwmethod", "nmulti", "ckertype", "np_tree", "seed_policy", "seed", "iter")
  keep[keep %in% names(a) & keep %in% names(b)]
}

summ_timing_by_combo <- function(df, label) {
  if (nrow(df) == 0L) return(data.frame())
  ks <- combo_keys(df)
  if (length(ks) == 0L) {
    return(data.frame(
      runs = nrow(df),
      label = label,
      mean_elapsed_bw = mean_or_na(df$elapsed_bw),
      median_elapsed_bw = median_or_na(df$elapsed_bw),
      mean_elapsed_fit = mean_or_na(df$elapsed_fit),
      median_elapsed_fit = median_or_na(df$elapsed_fit),
      mean_elapsed_total = mean_or_na(df$elapsed_total),
      median_elapsed_total = median_or_na(df$elapsed_total),
      stringsAsFactors = FALSE
    ))
  }
  split_ix <- interaction(df[ks], drop = TRUE, lex.order = TRUE)
  parts <- split(df, split_ix)
  rows <- lapply(parts, function(g) {
    data.frame(
      n = g$n[1],
      bwmethod = g$bwmethod[1],
      nmulti = if ("nmulti" %in% names(g)) g$nmulti[1] else NA_integer_,
      ckertype = g$ckertype[1],
      np_tree = g$np_tree[1],
      seed_policy = g$seed_policy[1],
      runs = nrow(g),
      label = label,
      mean_elapsed_bw = mean_or_na(g$elapsed_bw),
      median_elapsed_bw = median_or_na(g$elapsed_bw),
      mean_elapsed_fit = mean_or_na(g$elapsed_fit),
      median_elapsed_fit = median_or_na(g$elapsed_fit),
      mean_elapsed_total = mean_or_na(g$elapsed_total),
      median_elapsed_total = median_or_na(g$elapsed_total),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

mk_timing <- function(df, label) {
  data.frame(
    function_name = c("npplregbw", "npplreg", "npplreg_total"),
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

getcol <- function(df, nm, default = NA) {
  if (nm %in% names(df)) df[[nm]] else rep(default, nrow(df))
}

main <- function() {
  cfg <- parse_args(commandArgs(trailingOnly = TRUE))
  a <- prep_raw(read.csv(cfg$raw_a, stringsAsFactors = FALSE))
  b <- prep_raw(read.csv(cfg$raw_b, stringsAsFactors = FALSE))

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

  combo_a <- summ_timing_by_combo(a, cfg$label_a)
  combo_b <- summ_timing_by_combo(b, cfg$label_b)
  combo_keys_common <- intersect(
    c("n", "bwmethod", "nmulti", "ckertype", "np_tree", "seed_policy"),
    intersect(names(combo_a), names(combo_b))
  )
  combo_m <- merge(combo_a, combo_b, by = combo_keys_common, suffixes = c("_a", "_b"), all = TRUE)
  combo_timing <- data.frame(
    n = getcol(combo_m, "n"),
    bwmethod = getcol(combo_m, "bwmethod"),
    nmulti = getcol(combo_m, "nmulti"),
    ckertype = getcol(combo_m, "ckertype"),
    np_tree = getcol(combo_m, "np_tree"),
    seed_policy = getcol(combo_m, "seed_policy"),
    label_a = cfg$label_a,
    runs_a = combo_m$runs_a,
    mean_elapsed_bw_a = combo_m$mean_elapsed_bw_a,
    median_elapsed_bw_a = combo_m$median_elapsed_bw_a,
    mean_elapsed_fit_a = combo_m$mean_elapsed_fit_a,
    median_elapsed_fit_a = combo_m$median_elapsed_fit_a,
    mean_elapsed_total_a = combo_m$mean_elapsed_total_a,
    median_elapsed_total_a = combo_m$median_elapsed_total_a,
    label_b = cfg$label_b,
    runs_b = combo_m$runs_b,
    mean_elapsed_bw_b = combo_m$mean_elapsed_bw_b,
    median_elapsed_bw_b = combo_m$median_elapsed_bw_b,
    mean_elapsed_fit_b = combo_m$mean_elapsed_fit_b,
    median_elapsed_fit_b = combo_m$median_elapsed_fit_b,
    mean_elapsed_total_b = combo_m$mean_elapsed_total_b,
    median_elapsed_total_b = combo_m$median_elapsed_total_b,
    pct_faster_mean_bw_b_vs_a = mapply(pct_faster, combo_m$mean_elapsed_bw_a, combo_m$mean_elapsed_bw_b),
    pct_faster_median_bw_b_vs_a = mapply(pct_faster, combo_m$median_elapsed_bw_a, combo_m$median_elapsed_bw_b),
    pct_faster_mean_fit_b_vs_a = mapply(pct_faster, combo_m$mean_elapsed_fit_a, combo_m$mean_elapsed_fit_b),
    pct_faster_median_fit_b_vs_a = mapply(pct_faster, combo_m$median_elapsed_fit_a, combo_m$median_elapsed_fit_b),
    pct_faster_mean_total_b_vs_a = mapply(pct_faster, combo_m$mean_elapsed_total_a, combo_m$mean_elapsed_total_b),
    pct_faster_median_total_b_vs_a = mapply(pct_faster, combo_m$median_elapsed_total_a, combo_m$median_elapsed_total_b),
    stringsAsFactors = FALSE
  )
  write.csv(combo_timing, cfg$out_combo_timing, row.names = FALSE)

  keys <- match_keys(a, b)
  m <- merge(a, b, by = keys, suffixes = c("_a", "_b"), all = TRUE)

  absdiff <- function(x, y) ifelse(is.na(x) | is.na(y), NA_real_, abs(x - y))
  m$diff_fval <- absdiff(m$fval_a, m$fval_b)
  m$diff_ifval <- absdiff(m$ifval_a, m$ifval_b)
  m$diff_num_fval <- absdiff(m$num_fval_a, m$num_fval_b)
  m$diff_log_lik <- absdiff(m$log_likelihood_a, m$log_likelihood_b)
  m$same_bw <- ifelse(is.na(m$bw_a) | is.na(m$bw_b), NA, m$bw_a == m$bw_b)
  m$same_ok <- ifelse(is.na(m$ok_a) | is.na(m$ok_b), NA, m$ok_a == m$ok_b)

  objective <- data.frame(
    label_a = cfg$label_a,
    label_b = cfg$label_b,
    rows_compared = nrow(m),
    ok_match_rate = mean(m$same_ok, na.rm = TRUE),
    bw_match_rate = mean(m$same_bw, na.rm = TRUE),
    max_abs_diff_fval = max_or_na(m$diff_fval),
    max_abs_diff_ifval = max_or_na(m$diff_ifval),
    max_abs_diff_num_fval = max_or_na(m$diff_num_fval),
    max_abs_diff_log_likelihood = max_or_na(m$diff_log_lik),
    stringsAsFactors = FALSE
  )
  write.csv(objective, cfg$out_objective, row.names = FALSE)

  cks <- combo_keys(m)
  if (length(cks) == 0L) cks <- intersect(c("n", "bwmethod", "ckertype", "np_tree", "seed_policy"), names(m))
  split_ix <- interaction(m[cks], drop = TRUE, lex.order = TRUE)
  grp <- split(m, split_ix)
  combo_obj_rows <- lapply(grp, function(g) {
    data.frame(
      n = if ("n" %in% names(g)) g$n[1] else NA_integer_,
      bwmethod = if ("bwmethod" %in% names(g)) g$bwmethod[1] else "",
      nmulti = if ("nmulti" %in% names(g)) g$nmulti[1] else NA_integer_,
      ckertype = if ("ckertype" %in% names(g)) g$ckertype[1] else "",
      np_tree = if ("np_tree" %in% names(g)) g$np_tree[1] else NA,
      seed_policy = if ("seed_policy" %in% names(g)) g$seed_policy[1] else "",
      label_a = cfg$label_a,
      label_b = cfg$label_b,
      rows_compared = nrow(g),
      ok_match_rate = mean(g$same_ok, na.rm = TRUE),
      bw_match_rate = mean(g$same_bw, na.rm = TRUE),
      mean_abs_diff_fval = mean_or_na(g$diff_fval),
      max_abs_diff_fval = max_or_na(g$diff_fval),
      mean_abs_diff_ifval = mean_or_na(g$diff_ifval),
      max_abs_diff_ifval = max_or_na(g$diff_ifval),
      mean_abs_diff_num_fval = mean_or_na(g$diff_num_fval),
      max_abs_diff_num_fval = max_or_na(g$diff_num_fval),
      mean_abs_diff_log_likelihood = mean_or_na(g$diff_log_lik),
      max_abs_diff_log_likelihood = max_or_na(g$diff_log_lik),
      stringsAsFactors = FALSE
    )
  })
  combo_objective <- do.call(rbind, combo_obj_rows) %||% data.frame()
  write.csv(combo_objective, cfg$out_combo_objective, row.names = FALSE)

  cat("timing_csv=", cfg$out_timing, "\n", sep = "")
  cat("objective_csv=", cfg$out_objective, "\n", sep = "")
  cat("combo_timing_csv=", cfg$out_combo_timing, "\n", sep = "")
  cat("combo_objective_csv=", cfg$out_combo_objective, "\n", sep = "")
}

main()
