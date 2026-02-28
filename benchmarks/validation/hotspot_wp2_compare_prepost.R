args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: Rscript hotspot_wp2_compare_prepost.R <pre_raw.csv> <post_raw.csv> [out_dir]")
}

pre_path <- args[[1L]]
post_path <- args[[2L]]
out_dir <- if (length(args) >= 3L) args[[3L]] else tempfile("np_hotspot_prepost_")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pre <- read.csv(pre_path, stringsAsFactors = FALSE)
post <- read.csv(post_path, stringsAsFactors = FALSE)

key <- c("scenario_id", "function_name", "path", "regtype", "basis_degree", "n", "seed_policy", "iter", "seed")
missing_pre <- setdiff(key, names(pre))
missing_post <- setdiff(key, names(post))
if (length(missing_pre) || length(missing_post)) {
  stop(sprintf("Missing keys: pre=[%s], post=[%s]", paste(missing_pre, collapse = ","), paste(missing_post, collapse = ",")))
}

m <- merge(pre, post, by = key, suffixes = c("_pre", "_post"), all = FALSE)

is_true <- function(x) {
  if (is.logical(x)) return(isTRUE(x))
  if (is.numeric(x)) return(!is.na(x) && x != 0)
  if (is.character(x)) return(tolower(trimws(x)) %in% c("true", "t", "1", "yes"))
  FALSE
}

m$ok_pre_l <- vapply(m$ok_pre, is_true, logical(1))
m$ok_post_l <- vapply(m$ok_post, is_true, logical(1))

num_eq <- function(a, b, tol = 1e-12) {
  ok <- !is.na(a) & !is.na(b)
  out <- rep(FALSE, length(a))
  out[ok] <- abs(a[ok] - b[ok]) <= tol
  out
}

str_eq <- function(a, b) {
  ifelse(is.na(a) | is.na(b), FALSE, a == b)
}

m$num_feval_eq <- ifelse(is.na(m$num_feval_pre) & is.na(m$num_feval_post), TRUE,
                         ifelse(is.na(m$num_feval_pre) | is.na(m$num_feval_post), FALSE,
                                m$num_feval_pre == m$num_feval_post))
m$fval_eq <- ifelse(is.na(m$fval_pre) & is.na(m$fval_post), TRUE,
                    num_eq(m$fval_pre, m$fval_post, tol = 1e-12))
m$ifval_eq <- ifelse(is.na(m$ifval_pre) & is.na(m$ifval_post), TRUE,
                     ifelse(is.na(m$ifval_pre) | is.na(m$ifval_post), FALSE,
                            m$ifval_pre == m$ifval_post))
m$obj_hist_len_eq <- ifelse(is.na(m$obj_hist_len_pre) & is.na(m$obj_hist_len_post), TRUE,
                            ifelse(is.na(m$obj_hist_len_pre) | is.na(m$obj_hist_len_post), FALSE,
                                   m$obj_hist_len_pre == m$obj_hist_len_post))
m$obj_hist_sig_eq <- ifelse(is.na(m$obj_hist_sig_pre) & is.na(m$obj_hist_sig_post), TRUE,
                            str_eq(m$obj_hist_sig_pre, m$obj_hist_sig_post))
m$bw_sig_eq <- ifelse(is.na(m$bw_sig_pre) & is.na(m$bw_sig_post), TRUE,
                      str_eq(m$bw_sig_pre, m$bw_sig_post))
m$out_sig_eq <- ifelse(is.na(m$out_sig_pre) & is.na(m$out_sig_post), TRUE,
                       str_eq(m$out_sig_pre, m$out_sig_post))

m$elapsed_ratio_post_pre <- m$elapsed_sec_post / m$elapsed_sec_pre
m$elapsed_pct_change <- 100 * (m$elapsed_ratio_post_pre - 1)
m$gc_vcells_pct_change <- ifelse(is.na(m$gc_vcells_peak_mb_pre) | is.na(m$gc_vcells_peak_mb_post) | m$gc_vcells_peak_mb_pre == 0,
                                 NA_real_,
                                 100 * ((m$gc_vcells_peak_mb_post / m$gc_vcells_peak_mb_pre) - 1))

m$cv_path_lock <- ifelse(
  m$path == "bw",
  m$num_feval_eq & m$fval_eq & m$obj_hist_len_eq & m$obj_hist_sig_eq & m$bw_sig_eq,
  NA
)

m$parity_status <- ifelse(
  !(m$ok_pre_l & m$ok_post_l),
  "FAIL",
  ifelse(
    m$path == "bw" & !m$cv_path_lock,
    "PATH_DIVERGENCE",
    ifelse(m$bw_sig_eq & m$out_sig_eq & m$fval_eq, "PASS_STRICT", "FAIL")
  )
)

agg_keys <- c("scenario_id", "function_name", "path", "regtype", "basis_degree", "n", "seed_policy")

split_idx <- split(seq_len(nrow(m)), interaction(m[agg_keys], drop = TRUE))

agg_rows <- lapply(split_idx, function(idx) {
  d <- m[idx, , drop = FALSE]
  data.frame(
    scenario_id = d$scenario_id[1],
    function_name = d$function_name[1],
    path = d$path[1],
    regtype = d$regtype[1],
    basis_degree = d$basis_degree[1],
    n = d$n[1],
    seed_policy = d$seed_policy[1],
    runs = nrow(d),
    pre_mean_sec = mean(d$elapsed_sec_pre),
    post_mean_sec = mean(d$elapsed_sec_post),
    pre_median_sec = median(d$elapsed_sec_pre),
    post_median_sec = median(d$elapsed_sec_post),
    mean_ratio_post_pre = mean(d$elapsed_ratio_post_pre),
    mean_pct_change = 100 * (mean(d$elapsed_sec_post) / mean(d$elapsed_sec_pre) - 1),
    median_ratio_post_pre = median(d$elapsed_ratio_post_pre),
    median_pct_change = 100 * (median(d$elapsed_sec_post) / median(d$elapsed_sec_pre) - 1),
    mean_gc_vcells_pct_change = mean(d$gc_vcells_pct_change, na.rm = TRUE),
    parity_status = if (any(d$parity_status == "FAIL")) {
      "FAIL"
    } else if (any(d$parity_status == "PATH_DIVERGENCE")) {
      "PATH_DIVERGENCE"
    } else {
      "PASS_STRICT"
    },
    cv_path_lock_all = if (all(is.na(d$cv_path_lock))) NA else all(d$cv_path_lock, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
})

agg <- do.call(rbind, agg_rows)
agg <- agg[order(agg$function_name, agg$path, agg$regtype, agg$n, agg$seed_policy), ]

write.csv(m, file.path(out_dir, "wp2_prepost_rowlevel.csv"), row.names = FALSE)
write.csv(agg, file.path(out_dir, "wp2_prepost_summary.csv"), row.names = FALSE)

cat("pre:", pre_path, "\n")
cat("post:", post_path, "\n")
cat("out_dir:", out_dir, "\n")
cat("rowlevel:", file.path(out_dir, "wp2_prepost_rowlevel.csv"), "\n")
cat("summary:", file.path(out_dir, "wp2_prepost_summary.csv"), "\n")
cat("status counts:\n")
print(table(agg$parity_status, useNA = "ifany"))
cat("\nsummary preview:\n")
print(utils::head(agg, 16L))
