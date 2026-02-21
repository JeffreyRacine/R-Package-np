#!/usr/bin/env Rscript

parse_args <- function(args) {
  out <- list(
    timing_csv = "",
    objective_csv = "",
    combo_timing_csv = "",
    combo_objective_csv = "",
    out_md = "/tmp/npplreg_compare_report.md",
    title = "npplreg Comparison Report"
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
  req <- c("timing_csv", "objective_csv", "combo_timing_csv", "combo_objective_csv")
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

combo_cols <- function(df) {
  preferred <- c("regtype", "bwmethod", "ckertype", "np_tree", "seed_policy")
  keep <- preferred[preferred %in% names(df)]
  c(keep, "pct_faster_mean_total_b_vs_a", "pct_faster_median_total_b_vs_a")
}

main <- function() {
  cfg <- parse_args(commandArgs(trailingOnly = TRUE))
  timing <- read.csv(cfg$timing_csv, stringsAsFactors = FALSE)
  objective <- read.csv(cfg$objective_csv, stringsAsFactors = FALSE)
  combo_timing <- read.csv(cfg$combo_timing_csv, stringsAsFactors = FALSE)
  combo_objective <- read.csv(cfg$combo_objective_csv, stringsAsFactors = FALSE)

  ord <- order(combo_timing$pct_faster_mean_total_b_vs_a, decreasing = TRUE)
  ccols <- combo_cols(combo_timing)
  top_gain <- combo_timing[ord, ccols, drop = FALSE]
  top_loss <- combo_timing[rev(ord), ccols, drop = FALSE]
  top_gain <- head(top_gain, 5)
  top_loss <- head(top_loss, 5)

  timing2 <- timing
  for (nm in c("mean_a","median_a","mean_b","median_b","pct_faster_mean_b_vs_a","pct_faster_median_b_vs_a")) {
    timing2[[nm]] <- fmt_num(timing2[[nm]], d = 6L)
  }

  objective2 <- objective
  for (nm in setdiff(names(objective2), c("label_a","label_b"))) {
    objective2[[nm]] <- fmt_num(objective2[[nm]], d = 6L)
  }

  for (nm in c("pct_faster_mean_total_b_vs_a","pct_faster_median_total_b_vs_a")) {
    top_gain[[nm]] <- fmt_num(top_gain[[nm]], d = 6L)
    top_loss[[nm]] <- fmt_num(top_loss[[nm]], d = 6L)
  }

  lines <- c(
    paste0("# ", cfg$title),
    "",
    "## Inputs",
    paste0("- timing: `", cfg$timing_csv, "`"),
    paste0("- objective: `", cfg$objective_csv, "`"),
    paste0("- combo timing: `", cfg$combo_timing_csv, "`"),
    paste0("- combo objective: `", cfg$combo_objective_csv, "`"),
    "",
    "## Overall Timing",
    mk_table(timing2),
    "",
    "## Overall Objective Parity",
    mk_table(objective2),
    "",
    "## Top 5 Combo Gains (mean total)",
    mk_table(top_gain),
    "",
    "## Top 5 Combo Regressions (mean total)",
    mk_table(top_loss),
    "",
    "## Notes",
    "- Positive percent means `label_b` is faster than `label_a`.",
    "- Low bandwidth match rate can occur with optimizer path differences despite tight objective parity.",
    paste0("- Combo objective details are in `", cfg$combo_objective_csv, "`."),
    ""
  )

  writeLines(lines, con = cfg$out_md)
  cat("report_md=", cfg$out_md, "\n", sep = "")
}

main()
