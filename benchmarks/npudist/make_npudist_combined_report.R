#!/usr/bin/env Rscript

parse_args <- function(args) {
  out <- list(
    np_timing_csv = "",
    np_objective_csv = "",
    np_combo_timing_csv = "",
    np_combo_objective_csv = "",
    nprmpi_timing_csv = "",
    nprmpi_objective_csv = "",
    nprmpi_combo_timing_csv = "",
    nprmpi_combo_objective_csv = "",
    out_md = "/tmp/npudist_combined_compare_report.md",
    title = "npudist Combined Comparison Report"
  )

  for (a in args) {
    if (!grepl("^--", a)) stop("Bad arg: ", a)
    raw <- sub("^--", "", a)
    if (grepl("=", raw, fixed = TRUE)) {
      key <- sub("=.*$", "", raw)
      val <- sub("^[^=]*=", "", raw)
    } else {
      key <- raw
      val <- ""
    }
    if (key %in% names(out)) out[[key]] <- val else stop("Unknown arg: ", key)
  }

  req <- setdiff(names(out), c("out_md", "title"))
  miss <- req[!nzchar(unlist(out[req]))]
  if (length(miss)) stop("Missing required args: ", paste(miss, collapse = ", "))

  out
}

fmt_num <- function(x, d = 6L) {
  ifelse(is.na(x), "NA", formatC(x, digits = d, format = "fg", flag = "#"))
}

mk_table <- function(df) {
  cols <- names(df)
  header <- paste0("| ", paste(cols, collapse = " | "), " |")
  sep <- paste0("|", paste(rep("---", length(cols)), collapse = "|"), "|")
  rows <- apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  c(header, sep, rows)
}

top_combo <- function(df, n = 5L, best = TRUE) {
  ord <- order(df$pct_faster_mean_total_b_vs_a, decreasing = best)
  preferred <- c("regtype", "bwmethod", "ckertype", "np_tree", "seed_policy")
  keep <- preferred[preferred %in% names(df)]
  out <- df[ord, c(keep, "pct_faster_mean_total_b_vs_a", "pct_faster_median_total_b_vs_a"), drop = FALSE]
  head(out, n)
}

prep_timing <- function(df, backend_name) {
  keep <- c("function_name", "label_a", "mean_a", "median_a", "label_b", "mean_b", "median_b",
            "pct_faster_mean_b_vs_a", "pct_faster_median_b_vs_a")
  df <- df[, keep]
  df$backend <- backend_name
  df <- df[, c("backend", keep)]
  num_cols <- setdiff(names(df), c("backend", "function_name", "label_a", "label_b"))
  for (nm in num_cols) df[[nm]] <- fmt_num(df[[nm]])
  df
}

prep_objective <- function(df, backend_name) {
  df$backend <- backend_name
  df <- df[, c("backend", names(df)[names(df) != "backend"])]
  num_cols <- setdiff(names(df), c("backend", "label_a", "label_b"))
  for (nm in num_cols) df[[nm]] <- fmt_num(df[[nm]])
  df
}

main <- function() {
  cfg <- parse_args(commandArgs(trailingOnly = TRUE))

  np_timing <- read.csv(cfg$np_timing_csv, stringsAsFactors = FALSE)
  np_obj <- read.csv(cfg$np_objective_csv, stringsAsFactors = FALSE)
  np_combo_timing <- read.csv(cfg$np_combo_timing_csv, stringsAsFactors = FALSE)

  nprmpi_timing <- read.csv(cfg$nprmpi_timing_csv, stringsAsFactors = FALSE)
  nprmpi_obj <- read.csv(cfg$nprmpi_objective_csv, stringsAsFactors = FALSE)
  nprmpi_combo_timing <- read.csv(cfg$nprmpi_combo_timing_csv, stringsAsFactors = FALSE)

  timing_all <- rbind(prep_timing(np_timing, "np"), prep_timing(nprmpi_timing, "npRmpi"))
  objective_all <- rbind(prep_objective(np_obj, "np"), prep_objective(nprmpi_obj, "npRmpi"))

  np_gain <- top_combo(np_combo_timing, best = TRUE)
  np_loss <- top_combo(np_combo_timing, best = FALSE)
  nprmpi_gain <- top_combo(nprmpi_combo_timing, best = TRUE)
  nprmpi_loss <- top_combo(nprmpi_combo_timing, best = FALSE)

  for (nm in c("pct_faster_mean_total_b_vs_a", "pct_faster_median_total_b_vs_a")) {
    np_gain[[nm]] <- fmt_num(np_gain[[nm]])
    np_loss[[nm]] <- fmt_num(np_loss[[nm]])
    nprmpi_gain[[nm]] <- fmt_num(nprmpi_gain[[nm]])
    nprmpi_loss[[nm]] <- fmt_num(nprmpi_loss[[nm]])
  }

  lines <- c(
    paste0("# ", cfg$title),
    "",
    "## Inputs",
    paste0("- np timing: `", cfg$np_timing_csv, "`"),
    paste0("- np objective: `", cfg$np_objective_csv, "`"),
    paste0("- np combo timing: `", cfg$np_combo_timing_csv, "`"),
    paste0("- np combo objective: `", cfg$np_combo_objective_csv, "`"),
    paste0("- npRmpi timing: `", cfg$nprmpi_timing_csv, "`"),
    paste0("- npRmpi objective: `", cfg$nprmpi_objective_csv, "`"),
    paste0("- npRmpi combo timing: `", cfg$nprmpi_combo_timing_csv, "`"),
    paste0("- npRmpi combo objective: `", cfg$nprmpi_combo_objective_csv, "`"),
    "",
    "## Overall Timing (Both Backends)",
    mk_table(timing_all),
    "",
    "## Overall Objective Parity (Both Backends)",
    mk_table(objective_all),
    "",
    "## Top 5 Combo Gains: np",
    mk_table(np_gain),
    "",
    "## Top 5 Combo Regressions: np",
    mk_table(np_loss),
    "",
    "## Top 5 Combo Gains: npRmpi",
    mk_table(nprmpi_gain),
    "",
    "## Top 5 Combo Regressions: npRmpi",
    mk_table(nprmpi_loss),
    "",
    "## Notes",
    "- Positive percent means label_b is faster than label_a.",
    "- Bandwidth match rate can be low when optimizer paths differ despite small objective differences.",
    "- Use backend-specific combo objective CSV files for per-combo parity details.",
    ""
  )

  writeLines(lines, con = cfg$out_md)
  cat("combined_report_md=", cfg$out_md, "\n", sep = "")
}

main()
