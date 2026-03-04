#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  i <- match(name, args)
  if (is.na(i) || i == length(args)) return(default)
  args[[i + 1L]]
}

repo <- get_arg("--repo", "/Users/jracine/Development/np-master")
n <- as.integer(get_arg("--n", "1000"))
seed <- as.integer(get_arg("--seed", "42"))
boot_num <- as.integer(get_arg("--boot-num", "399"))
casefile <- get_arg("--case-file", "")
out_dir <- get_arg("--out-dir", ".")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(devtools))
load_all(repo, quiet = TRUE, export_all = FALSE)

build_default_case <- function(n, seed) {
  set.seed(seed)
  x1 <- runif(n)
  x2 <- runif(n)
  f1 <- factor(sample(c("a", "b", "c", "d"), n, replace = TRUE))
  xdat <- data.frame(x1 = x1, x2 = x2, f1 = f1)
  ydat <- sin(2 * pi * x1) + 0.5 * x2 + as.numeric(f1) / 10 + rnorm(n, sd = 0.2)
  bws <- npregbw(xdat = xdat, ydat = ydat, regtype = "lc", nmulti = 1)
  list(xdat = xdat, ydat = ydat, bws = bws)
}

if (nzchar(casefile)) {
  env <- new.env(parent = .GlobalEnv)
  sys.source(casefile, envir = env)
  if (!exists("build_case", envir = env, inherits = FALSE))
    stop("case file must define build_case(n, seed)")
  cs <- env$build_case(n = n, seed = seed)
} else {
  cs <- build_default_case(n = n, seed = seed)
}

xdat <- cs$xdat
ydat <- cs$ydat
bws <- cs$bws

has_nan_inf_core <- function(x) {
  if (is.numeric(x))
    return(any(is.nan(x) | is.infinite(x)))

  if (is.list(x)) {
    nms <- names(x)
    for (i in seq_along(x)) {
      nm <- if (!is.null(nms) && nzchar(nms[i])) nms[i] else ""
      # Skip embedded bandwidth metadata (contains intentional +/-Inf kernel bounds).
      if (identical(nm, "bws"))
        next
      if (has_nan_inf_core(x[[i]]))
        return(TRUE)
    }
  }

  FALSE
}

status <- "ok"
err <- ""
nan_inf <- FALSE
elapsed <- NA_real_
warn_log <- character()

start_time <- proc.time()[["elapsed"]]
res <- tryCatch(
  withCallingHandlers(
    plot(
      bws,
      xdat = xdat,
      ydat = ydat,
      plot.behavior = "data",
      gradients = FALSE,
      neval = 100,
      plot.errors.method = "bootstrap",
      plot.errors.boot.method = "wild",
      plot.errors.boot.wild = "rademacher",
      plot.errors.boot.num = boot_num,
      plot.errors.type = "pointwise",
      plot.errors.center = "estimate"
    ),
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
elapsed <- proc.time()[["elapsed"]] - start_time

if (identical(status, "ok") && is.null(res))
  status <- "null_return"

if (!is.null(res))
  nan_inf <- isTRUE(has_nan_inf_core(res))

warn_file <- file.path(out_dir, "warnings.log")
if (length(warn_log)) writeLines(warn_log, warn_file)

cat(sprintf("status=%s\n", status))
cat(sprintf("n=%d\n", n))
cat(sprintf("seed=%d\n", seed))
cat(sprintf("boot_num=%d\n", boot_num))
cat(sprintf("elapsed_sec=%.6f\n", elapsed))
cat(sprintf("nan_inf=%s\n", if (isTRUE(nan_inf)) "TRUE" else "FALSE"))
cat(sprintf("warn_count=%d\n", length(warn_log)))
if (nzchar(err)) cat(sprintf("error=%s\n", err))

if (!identical(status, "ok") || isTRUE(nan_inf))
  quit(save = "no", status = 10)

quit(save = "no", status = 0)
