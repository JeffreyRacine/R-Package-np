#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  i <- match(name, args)
  if (is.na(i) || i == length(args)) return(default)
  args[[i + 1L]]
}

repo <- get_arg("--repo", "/Users/jracine/Development/np-master")
label <- get_arg("--label", "candidate")
out_csv <- get_arg("--out-csv", tempfile("goal1_perf_runs_", fileext = ".csv"))
boot_num <- as.integer(get_arg("--boot-num", "79"))
blocklen <- as.integer(get_arg("--blocklen", "5"))
seed_text <- get_arg("--seeds", "")

if (!nzchar(seed_text)) {
  seeds <- 1:25
} else {
  seeds <- as.integer(strsplit(seed_text, ",", fixed = TRUE)[[1L]])
  seeds <- seeds[is.finite(seeds)]
  if (!length(seeds)) stop("no valid seeds provided")
}

dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)

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

set.seed(20260304)
n <- 150L
x <- runif(n)
x2 <- runif(n)
f <- factor(sample(c("a", "b", "c", "d"), n, replace = TRUE))
y <- sin(2 * pi * x) + 0.4 * x2 + as.numeric(f) / 9 + rnorm(n, sd = 0.12)

x_reg <- data.frame(x = x, f = f)
reg_bw <- npregbw(xdat = x_reg, ydat = y, regtype = "lc", bws = c(0.35, 0.5), bandwidth.compute = FALSE)
u_bw <- npudensbw(dat = data.frame(x = x), bws = 0.35, bandwidth.compute = FALSE)

run_plot <- function(expr) {
  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
}

scenarios <- list(
  list(
    name = "regression_wild",
    fn = function() run_plot(plot(
      reg_bw, xdat = x_reg, ydat = y, plot.behavior = "data", perspective = FALSE,
      gradients = FALSE, plot.errors.method = "bootstrap", plot.errors.boot.method = "wild",
      plot.errors.boot.wild = "rademacher", plot.errors.boot.num = boot_num, plot.errors.boot.blocklen = blocklen
    ))
  ),
  list(
    name = "regression_inid",
    fn = function() run_plot(plot(
      reg_bw, xdat = x_reg, ydat = y, plot.behavior = "data", perspective = FALSE,
      gradients = FALSE, plot.errors.method = "bootstrap", plot.errors.boot.method = "inid",
      plot.errors.boot.num = boot_num, plot.errors.boot.blocklen = blocklen
    ))
  ),
  list(
    name = "regression_fixed",
    fn = function() run_plot(plot(
      reg_bw, xdat = x_reg, ydat = y, plot.behavior = "data", perspective = FALSE,
      gradients = FALSE, plot.errors.method = "bootstrap", plot.errors.boot.method = "fixed",
      plot.errors.boot.num = boot_num, plot.errors.boot.blocklen = blocklen
    ))
  ),
  list(
    name = "regression_geom",
    fn = function() run_plot(plot(
      reg_bw, xdat = x_reg, ydat = y, plot.behavior = "data", perspective = FALSE,
      gradients = FALSE, plot.errors.method = "bootstrap", plot.errors.boot.method = "geom",
      plot.errors.boot.num = boot_num, plot.errors.boot.blocklen = blocklen
    ))
  ),
  list(
    name = "udens_inid",
    fn = function() run_plot(plot(
      u_bw, plot.behavior = "data", perspective = FALSE,
      plot.errors.method = "bootstrap", plot.errors.boot.method = "inid",
      plot.errors.boot.num = boot_num, plot.errors.boot.blocklen = blocklen
    ))
  ),
  list(
    name = "udens_fixed",
    fn = function() run_plot(plot(
      u_bw, plot.behavior = "data", perspective = FALSE,
      plot.errors.method = "bootstrap", plot.errors.boot.method = "fixed",
      plot.errors.boot.num = boot_num, plot.errors.boot.blocklen = blocklen
    ))
  )
)

rows <- list()
k <- 1L
for (sd in seeds) {
  set.seed(sd)
  for (sc in scenarios) {
    warn_log <- character()
    status <- "ok"
    err <- ""
    out <- NULL
    elapsed <- system.time({
      out <- tryCatch(
        withCallingHandlers(
          sc$fn(),
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
    })[["elapsed"]]
    nonfinite <- if (!is.null(out)) isTRUE(has_nan_inf(out)) else FALSE

    rows[[k]] <- data.frame(
      label = label,
      scenario = sc$name,
      seed = as.integer(sd),
      elapsed_sec = as.numeric(elapsed),
      status = status,
      warn_count = as.integer(length(warn_log)),
      nonfinite = as.logical(nonfinite),
      error = err,
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
}

out <- do.call(rbind, rows)
write.csv(out, out_csv, row.names = FALSE)

bad <- subset(out, status != "ok" | nonfinite)
cat(sprintf("label=%s\n", label))
cat(sprintf("repo=%s\n", repo))
cat(sprintf("rows=%d\n", nrow(out)))
cat(sprintf("bad_rows=%d\n", nrow(bad)))
cat(sprintf("out_csv=%s\n", out_csv))

if (nrow(bad) > 0L) quit(save = "no", status = 1L)
quit(save = "no", status = 0L)
