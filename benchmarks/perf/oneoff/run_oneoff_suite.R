#!/usr/bin/env Rscript

parse_args <- function(args) {
  out <- list(
    n_values = c(100L, 250L, 500L),
    times = 50L,
    seed_policy = "varying",
    base_seed = 42L,
    out_manifest = "/tmp/np_oneoff_suite_manifest.csv",
    show_progress = TRUE
  )
  if (length(args) == 0L) return(out)
  for (a in args) {
    if (!grepl("^--", a)) stop("Bad arg: ", a)
    raw <- sub("^--", "", a)
    key <- sub("=.*$", "", raw)
    val <- if (grepl("=", raw, fixed = TRUE)) sub("^[^=]*=", "", raw) else "TRUE"
    if (key == "n_values") out$n_values <- as.integer(strsplit(val, ",", fixed = TRUE)[[1]])
    else if (key == "times") out$times <- as.integer(val)
    else if (key == "seed_policy") out$seed_policy <- val
    else if (key == "base_seed") out$base_seed <- as.integer(val)
    else if (key == "out_manifest") out$out_manifest <- val
    else if (key == "show_progress") out$show_progress <- as.logical(val)
    else stop("Unknown arg: ", key)
  }
  out
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  cfg <- parse_args(args)
  script <- "/Users/jracine/Development/np-master/benchmarks/perf/oneoff/bench_oneoff_param.R"

  real_funs <- c("npcmstest", "npconmode", "npqreg")
  synth_funs <- c("npcopula", "npdeneqtest", "npdeptest", "npregiv", "npsdeptest", "npsigtest", "npsymtest", "npunitest")

  rows <- list()
  idx <- 0L

  for (fun in real_funs) {
    idx <- idx + 1L
    out_raw <- sprintf("/tmp/np_oneoff_%s_raw.csv", fun)
    out_sum <- sprintf("/tmp/np_oneoff_%s_summary.csv", fun)
    cmd <- c(script,
             paste0("--fun=", fun),
             "--n=100",
             paste0("--times=", cfg$times),
             paste0("--seed_policy=", cfg$seed_policy),
             paste0("--base_seed=", cfg$base_seed),
             paste0("--out_raw=", out_raw),
             paste0("--out_summary=", out_sum),
             "--show_progress=FALSE")
    if (cfg$show_progress) cat("running", fun, "(real data)\n")
    rc <- system2("Rscript", cmd)
    rows[[idx]] <- data.frame(fun = fun, n = NA_integer_, out_raw = out_raw, out_summary = out_sum, rc = rc, stringsAsFactors = FALSE)
  }

  for (fun in synth_funs) {
    for (n in cfg$n_values) {
      idx <- idx + 1L
      out_raw <- sprintf("/tmp/np_oneoff_%s_n%d_raw.csv", fun, n)
      out_sum <- sprintf("/tmp/np_oneoff_%s_n%d_summary.csv", fun, n)
      cmd <- c(script,
               paste0("--fun=", fun),
               paste0("--n=", n),
               paste0("--times=", cfg$times),
               paste0("--seed_policy=", cfg$seed_policy),
               paste0("--base_seed=", cfg$base_seed),
               paste0("--out_raw=", out_raw),
               paste0("--out_summary=", out_sum),
               "--show_progress=FALSE")
      if (cfg$show_progress) cat("running", fun, "n=", n, "\n")
      rc <- system2("Rscript", cmd)
      rows[[idx]] <- data.frame(fun = fun, n = n, out_raw = out_raw, out_summary = out_sum, rc = rc, stringsAsFactors = FALSE)
    }
  }

  m <- do.call(rbind, rows)
  write.csv(m, cfg$out_manifest, row.names = FALSE)
  if (cfg$show_progress) {
    cat("manifest:", cfg$out_manifest, "\n")
    print(m)
  }
}

if (sys.nframe() == 0) main()
