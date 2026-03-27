#!/usr/bin/env Rscript

parse_args <- function(args) {
  out <- list(
    n_values = c(100L, 250L, 500L),
    nslaves = 1L,
    times = 50L,
    seed_policy = "varying",
    base_seed = 42L,
    iface = "en0",
    out_manifest = "/tmp/nprmpi_oneoff_suite_manifest.csv",
    show_progress = TRUE
  )
  if (length(args) == 0L) return(out)
  for (a in args) {
    if (!grepl("^--", a)) stop("Bad arg: ", a)
    raw <- sub("^--", "", a)
    key <- sub("=.*$", "", raw)
    val <- if (grepl("=", raw, fixed = TRUE)) sub("^[^=]*=", "", raw) else "TRUE"
    if (key == "n_values") out$n_values <- as.integer(strsplit(val, ",", fixed = TRUE)[[1]])
    else if (key == "nslaves" || key == "rslaves") out$nslaves <- as.integer(val)
    else if (key == "times") out$times <- as.integer(val)
    else if (key == "seed_policy") out$seed_policy <- val
    else if (key == "base_seed") out$base_seed <- as.integer(val)
    else if (key == "iface") out$iface <- val
    else if (key == "out_manifest") out$out_manifest <- val
    else if (key == "show_progress") out$show_progress <- as.logical(val)
    else stop("Unknown arg: ", key)
  }
  out
}

blas_string <- function() {
  blas <- sessionInfo()$BLAS
  if (is.null(blas) || length(blas) == 0L || is.na(blas[1L]) || !nzchar(blas[1L])) {
    return("unknown")
  }
  as.character(blas[1L])
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  cfg <- parse_args(args)
  script <- "/Users/jracine/Development/np-npRmpi/benchmarks/perf/oneoff/bench_oneoff_param_nprmpi.R"
  log_dir <- dirname(cfg$out_manifest)
  provenance_file <- paste0(sub("\\.csv$", "", cfg$out_manifest), "_provenance.txt")

  real_funs <- c("npcmstest", "npconmode", "npqreg")
  synth_funs <- c("npcopula", "npdeneqtest", "npdeptest", "npregiv", "npsdeptest", "npsigtest", "npsymtest", "npunitest")

  rows <- list()
  idx <- 0L

  run_case <- function(fun, n_value, is_real) {
    out_raw <- if (is_real) sprintf("/tmp/nprmpi_oneoff_%s_raw.csv", fun) else sprintf("/tmp/nprmpi_oneoff_%s_n%d_raw.csv", fun, n_value)
    out_sum <- if (is_real) sprintf("/tmp/nprmpi_oneoff_%s_summary.csv", fun) else sprintf("/tmp/nprmpi_oneoff_%s_n%d_summary.csv", fun, n_value)
    log_file <- file.path(log_dir, if (is_real) sprintf("oneoff_%s.log", fun) else sprintf("oneoff_%s_n%d.log", fun, n_value))

    args <- c(script,
              paste0("--fun=", fun),
              paste0("--n=", if (is_real) 100L else n_value),
              paste0("--nslaves=", cfg$nslaves),
              paste0("--times=", cfg$times),
              paste0("--seed_policy=", cfg$seed_policy),
              paste0("--base_seed=", cfg$base_seed),
              paste0("--out_raw=", out_raw),
              paste0("--out_summary=", out_sum),
              "--show_progress=FALSE")

    env <- c(sprintf("FI_TCP_IFACE=%s", cfg$iface))
    rc <- system2("env", c(env, "script", "-q", "/dev/null", "Rscript", args), stdout = log_file, stderr = log_file)
    data.frame(fun = fun, n = if (is_real) NA_integer_ else n_value, out_raw = out_raw, out_summary = out_sum, log = log_file, rc = rc, stringsAsFactors = FALSE)
  }

  for (fun in real_funs) {
    idx <- idx + 1L
    if (cfg$show_progress) cat("running", fun, "(real data)\n")
    rows[[idx]] <- run_case(fun, NA_integer_, TRUE)
  }

  for (fun in synth_funs) {
    for (n in cfg$n_values) {
      idx <- idx + 1L
      if (cfg$show_progress) cat("running", fun, "n=", n, "\n")
      rows[[idx]] <- run_case(fun, n, FALSE)
    }
  }

  m <- do.call(rbind, rows)
  writeLines(
    c(
      paste0("R.version.string=", R.version.string),
      paste0("BLAS=", blas_string())
    ),
    con = provenance_file
  )
  write.csv(m, cfg$out_manifest, row.names = FALSE)
  if (cfg$show_progress) {
    cat("R.version.string=", R.version.string, "\n", sep = "")
    cat("BLAS=", blas_string(), "\n", sep = "")
    cat("provenance:", provenance_file, "\n")
    cat("manifest:", cfg$out_manifest, "\n")
    print(m)
  }
}

if (sys.nframe() == 0) main()
