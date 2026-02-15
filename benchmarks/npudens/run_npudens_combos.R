#!/usr/bin/env Rscript

parse_args <- function(args) {
  out <- list(
    n = 100L,
    times = 50L,
    base_seed = 42L,
    nmulti = 1L,
    out_dir = "/tmp",
    tag = "",
    show_progress = TRUE
  )

  if (length(args) == 0L) return(out)

  for (a in args) {
    if (!grepl("^--", a)) stop("Bad arg: ", a)
    raw <- sub("^--", "", a)
    key <- sub("=.*$", "", raw)
    val <- if (grepl("=", raw, fixed = TRUE)) sub("^[^=]*=", "", raw) else "TRUE"

    if (key == "n") out$n <- as.integer(val)
    else if (key == "times") out$times <- as.integer(val)
    else if (key == "base_seed") out$base_seed <- as.integer(val)
    else if (key == "nmulti") out$nmulti <- as.integer(val)
    else if (key == "out_dir") out$out_dir <- val
    else if (key == "tag") out$tag <- val
    else if (key == "show_progress") out$show_progress <- as.logical(val)
    else stop("Unknown arg: ", key)
  }

  out
}

repo_hash <- function(repo_dir) {
  x <- tryCatch(system2("git", c("-C", repo_dir, "rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE), error = function(e) "")
  x <- trimws(paste(x, collapse = ""))
  if (nzchar(x)) x else "nogit"
}

mk_run_id <- function(cfg, backend, hash) {
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  tag <- if (nzchar(cfg$tag)) paste0("_", cfg$tag) else ""
  sprintf("%s_%s_%s_n%d_t%d_s%d_m%d%s", backend, hash, ts, cfg$n, cfg$times, cfg$base_seed, cfg$nmulti, tag)
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  cfg <- parse_args(args)

  argv <- commandArgs(trailingOnly = FALSE)
  script_arg <- argv[grep("^--file=", argv)][1]
  script_path <- if (!is.na(script_arg) && nzchar(script_arg)) sub("^--file=", "", script_arg) else file.path(getwd(), "run_npudens_combos.R")
  here <- dirname(normalizePath(script_path))

  bench_script <- file.path(here, "bench_npudens_param.R")
  if (!file.exists(bench_script)) stop("Missing benchmark script: ", bench_script)

  hash <- repo_hash(here)
  run_id <- mk_run_id(cfg, "np", hash)

  out_raw <- file.path(cfg$out_dir, paste0("npudens_combo_raw_", run_id, ".csv"))
  out_summary <- file.path(cfg$out_dir, paste0("npudens_combo_summary_", run_id, ".csv"))
  out_manifest <- file.path(cfg$out_dir, paste0("npudens_combo_manifest_", run_id, ".csv"))

  combos <- expand.grid(
    bwmethod = c("cv.ml", "cv.ls"),
    ckertype = c("gaussian", "epanechnikov"),
    np_tree = c(FALSE, TRUE),
    seed_policy = c("fixed", "varying"),
    stringsAsFactors = FALSE
  )

  status <- vector("list", nrow(combos))

  for (i in seq_len(nrow(combos))) {
    row <- combos[i, ]

    if (cfg$show_progress) {
      cat(sprintf("[%02d/%02d] bw=%s cker=%s tree=%s seed_policy=%s\n",
                  i, nrow(combos), row$bwmethod, row$ckertype,
                  as.character(row$np_tree), row$seed_policy))
    }

    cmd_args <- c(
      bench_script,
      sprintf("--n=%d", cfg$n),
      sprintf("--times=%d", cfg$times),
      sprintf("--base_seed=%d", cfg$base_seed),
      sprintf("--nmulti=%d", cfg$nmulti),
      sprintf("--bwmethod=%s", row$bwmethod),
      sprintf("--ckertype=%s", row$ckertype),
      sprintf("--np_tree=%s", toupper(as.character(row$np_tree))),
      sprintf("--seed_policy=%s", row$seed_policy),
      sprintf("--out_raw=%s", out_raw),
      sprintf("--out_summary=%s", out_summary),
      "--show_progress=FALSE"
    )

    rc <- system2("Rscript", cmd_args)

    status[[i]] <- data.frame(
      combo_id = i,
      bwmethod = row$bwmethod,
      ckertype = row$ckertype,
      np_tree = row$np_tree,
      seed_policy = row$seed_policy,
      rc = rc,
      stringsAsFactors = FALSE
    )
  }

  manifest <- do.call(rbind, status)
  write.csv(manifest, out_manifest, row.names = FALSE)

  cat("raw_csv=", out_raw, "\n", sep = "")
  cat("summary_csv=", out_summary, "\n", sep = "")
  cat("manifest_csv=", out_manifest, "\n", sep = "")

  if (any(manifest$rc != 0L)) stop("One or more combos failed. See manifest: ", out_manifest)
}

if (sys.nframe() == 0) main()
