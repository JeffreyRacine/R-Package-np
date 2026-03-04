#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  i <- match(name, args)
  if (is.na(i) || i == length(args)) return(default)
  args[[i + 1L]]
}

repo <- get_arg("--repo", "/Users/jracine/Development/np-master")
out_dir <- get_arg("--out-dir", tempfile("goal1_matrix_"))
boot_num <- as.integer(get_arg("--boot-num", "79"))
blocklen <- as.integer(get_arg("--blocklen", "5"))

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(devtools))
load_all(repo, quiet = TRUE, export_all = FALSE)

has_nan_inf <- function(x) {
  if (is.numeric(x)) return(any(is.nan(x) | is.infinite(x)))
  if (is.list(x)) {
    nms <- names(x)
    for (i in seq_along(x)) {
      nm <- if (!is.null(nms) && nzchar(nms[i])) nms[i] else ""
      if (identical(nm, "bws")) next
      if (has_nan_inf(x[[i]])) return(TRUE)
    }
  }
  FALSE
}

run_case <- function(case_name, expected, error_pattern = "", runner) {
  warn_log <- character()
  status <- "ok"
  err <- ""
  out <- NULL
  t0 <- proc.time()[["elapsed"]]

  out <- tryCatch(
    withCallingHandlers(
      runner(),
      warning = function(w) {
        warn_log <<- c(warn_log, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      status <<- "error"
      err <<- conditionMessage(e)
      NULL
    }
  )

  elapsed <- proc.time()[["elapsed"]] - t0
  nonfinite <- if (!is.null(out)) isTRUE(has_nan_inf(out)) else FALSE

  pass <- FALSE
  if (identical(expected, "ok")) {
    pass <- identical(status, "ok") && !isTRUE(nonfinite)
  } else if (identical(expected, "error")) {
    pass <- identical(status, "error") &&
      (if (nzchar(error_pattern)) isTRUE(grepl(error_pattern, err, fixed = TRUE)) else TRUE)
  } else {
    stop("invalid expected value")
  }

  list(
    scenario = case_name,
    expected = expected,
    status = status,
    pass = pass,
    warn_count = length(warn_log),
    nonfinite = nonfinite,
    elapsed_sec = elapsed,
    error = err,
    warnings = warn_log
  )
}

set.seed(20260304)
n <- 140L
x <- runif(n)
x2 <- runif(n)
f <- factor(sample(c("a", "b", "c", "d"), n, replace = TRUE))
y <- sin(2 * pi * x) + 0.4 * x2 + as.numeric(f) / 9 + rnorm(n, sd = 0.12)

x_reg <- data.frame(x = x, f = f)
reg_bw <- npregbw(xdat = x_reg, ydat = y, regtype = "lc", bws = c(0.35, 0.5), bandwidth.compute = FALSE)

u_bw <- npudensbw(dat = data.frame(x = x), bws = 0.35, bandwidth.compute = FALSE)

x_cd <- data.frame(x = x)
y_cd <- data.frame(y = x2)
c_bw <- npcdensbw(xdat = x_cd, ydat = y_cd, bws = c(0.35, 0.35), bandwidth.compute = FALSE)

run_plot <- function(expr) {
  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
}

cases <- list(
  list(
    name = "regression_wild",
    expected = "ok",
    pattern = "",
    fn = function() run_plot(plot(
      reg_bw, xdat = x_reg, ydat = y, plot.behavior = "data", perspective = FALSE,
      gradients = FALSE, plot.errors.method = "bootstrap", plot.errors.boot.method = "wild",
      plot.errors.boot.wild = "rademacher", plot.errors.boot.num = boot_num,
      plot.errors.boot.blocklen = blocklen
    ))
  ),
  list(
    name = "regression_inid",
    expected = "ok",
    pattern = "",
    fn = function() run_plot(plot(
      reg_bw, xdat = x_reg, ydat = y, plot.behavior = "data", perspective = FALSE,
      gradients = FALSE, plot.errors.method = "bootstrap", plot.errors.boot.method = "inid",
      plot.errors.boot.num = boot_num, plot.errors.boot.blocklen = blocklen
    ))
  ),
  list(
    name = "regression_fixed",
    expected = "ok",
    pattern = "",
    fn = function() run_plot(plot(
      reg_bw, xdat = x_reg, ydat = y, plot.behavior = "data", perspective = FALSE,
      gradients = FALSE, plot.errors.method = "bootstrap", plot.errors.boot.method = "fixed",
      plot.errors.boot.num = boot_num, plot.errors.boot.blocklen = blocklen
    ))
  ),
  list(
    name = "regression_geom",
    expected = "ok",
    pattern = "",
    fn = function() run_plot(plot(
      reg_bw, xdat = x_reg, ydat = y, plot.behavior = "data", perspective = FALSE,
      gradients = FALSE, plot.errors.method = "bootstrap", plot.errors.boot.method = "geom",
      plot.errors.boot.num = boot_num, plot.errors.boot.blocklen = blocklen
    ))
  ),
  list(
    name = "udens_wild_reject",
    expected = "error",
    pattern = "not supported for unconditional density/distribution estimators",
    fn = function() run_plot(plot(
      u_bw, plot.behavior = "data", perspective = FALSE,
      plot.errors.method = "bootstrap", plot.errors.boot.method = "wild",
      plot.errors.boot.num = boot_num, plot.errors.boot.blocklen = blocklen
    ))
  ),
  list(
    name = "udens_inid",
    expected = "ok",
    pattern = "",
    fn = function() run_plot(plot(
      u_bw, plot.behavior = "data", perspective = FALSE,
      plot.errors.method = "bootstrap", plot.errors.boot.method = "inid",
      plot.errors.boot.num = boot_num, plot.errors.boot.blocklen = blocklen
    ))
  ),
  list(
    name = "udens_fixed",
    expected = "ok",
    pattern = "",
    fn = function() run_plot(plot(
      u_bw, plot.behavior = "data", perspective = FALSE,
      plot.errors.method = "bootstrap", plot.errors.boot.method = "fixed",
      plot.errors.boot.num = boot_num, plot.errors.boot.blocklen = blocklen
    ))
  ),
  list(
    name = "udens_geom",
    expected = "ok",
    pattern = "",
    fn = function() run_plot(plot(
      u_bw, plot.behavior = "data", perspective = FALSE,
      plot.errors.method = "bootstrap", plot.errors.boot.method = "geom",
      plot.errors.boot.num = boot_num, plot.errors.boot.blocklen = blocklen
    ))
  ),
  list(
    name = "cdens_wild_reject",
    expected = "error",
    pattern = "not supported for conditional density/distribution estimators",
    fn = function() run_plot(plot(
      c_bw, xdat = x_cd, ydat = y_cd, plot.behavior = "data", perspective = FALSE,
      gradients = FALSE, plot.errors.method = "bootstrap", plot.errors.boot.method = "wild",
      plot.errors.boot.num = boot_num, plot.errors.boot.blocklen = blocklen
    ))
  ),
  list(
    name = "cdens_inid",
    expected = "ok",
    pattern = "",
    fn = function() run_plot(plot(
      c_bw, xdat = x_cd, ydat = y_cd, plot.behavior = "data", perspective = FALSE,
      gradients = FALSE, plot.errors.method = "bootstrap", plot.errors.boot.method = "inid",
      plot.errors.boot.num = boot_num, plot.errors.boot.blocklen = blocklen
    ))
  ),
  list(
    name = "cdens_fixed",
    expected = "ok",
    pattern = "",
    fn = function() run_plot(plot(
      c_bw, xdat = x_cd, ydat = y_cd, plot.behavior = "data", perspective = FALSE,
      gradients = FALSE, plot.errors.method = "bootstrap", plot.errors.boot.method = "fixed",
      plot.errors.boot.num = boot_num, plot.errors.boot.blocklen = blocklen
    ))
  ),
  list(
    name = "cdens_geom",
    expected = "ok",
    pattern = "",
    fn = function() run_plot(plot(
      c_bw, xdat = x_cd, ydat = y_cd, plot.behavior = "data", perspective = FALSE,
      gradients = FALSE, plot.errors.method = "bootstrap", plot.errors.boot.method = "geom",
      plot.errors.boot.num = boot_num, plot.errors.boot.blocklen = blocklen
    ))
  )
)

rows <- vector("list", length(cases))
warn_dir <- file.path(out_dir, "warnings")
dir.create(warn_dir, recursive = TRUE, showWarnings = FALSE)

for (i in seq_along(cases)) {
  cs <- cases[[i]]
  res <- run_case(
    case_name = cs$name,
    expected = cs$expected,
    error_pattern = cs$pattern,
    runner = cs$fn
  )
  rows[[i]] <- data.frame(
    scenario = res$scenario,
    expected = res$expected,
    status = res$status,
    pass = as.logical(res$pass),
    warn_count = as.integer(res$warn_count),
    nonfinite = as.logical(res$nonfinite),
    elapsed_sec = as.numeric(res$elapsed_sec),
    error = as.character(res$error),
    stringsAsFactors = FALSE
  )
  if (length(res$warnings)) {
    writeLines(res$warnings, con = file.path(warn_dir, paste0(cs$name, ".warnings.log")))
  }
}

out <- do.call(rbind, rows)
csv_path <- file.path(out_dir, "matrix_results.csv")
write.csv(out, csv_path, row.names = FALSE)

pass_all <- isTRUE(all(out$pass))
summary_path <- file.path(out_dir, "matrix_summary.txt")
writeLines(
  c(
    sprintf("repo=%s", repo),
    sprintf("boot_num=%d", boot_num),
    sprintf("blocklen=%d", blocklen),
    sprintf("cases=%d", nrow(out)),
    sprintf("pass_cases=%d", sum(out$pass)),
    sprintf("fail_cases=%d", sum(!out$pass)),
    sprintf("result=%s", if (pass_all) "PASS" else "FAIL"),
    sprintf("csv=%s", csv_path)
  ),
  con = summary_path
)

cat(sprintf("matrix_csv=%s\n", csv_path))
cat(sprintf("matrix_summary=%s\n", summary_path))
cat(sprintf("matrix_result=%s\n", if (pass_all) "PASS" else "FAIL"))

if (!pass_all) quit(save = "no", status = 1L)
quit(save = "no", status = 0L)
