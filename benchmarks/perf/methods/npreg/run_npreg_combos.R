#!/usr/bin/env Rscript

parse_args <- function(args) {
  parse_csv <- function(val) {
    trimws(strsplit(val, ",", fixed = TRUE)[[1]])
  }
  parse_csv_bool <- function(val) {
    v <- tolower(parse_csv(val))
    if (!all(v %in% c("true", "false"))) stop("boolean CSV fields must use TRUE/FALSE")
    v == "true"
  }
  out <- list(
    n = 100L,
    times = 50L,
    base_seed = 42L,
    nmulti = 1L,
    lp_bases = c("glp", "additive", "tensor"),
    lp_degree = "2,2",
    lp_bernstein.basis = c(FALSE, TRUE),
    nslaves = 1L,
    out_dir = "/tmp",
    fi_tcp_iface = "en0",
    r_libs = Sys.getenv("R_LIBS", unset = ""),
    max_combos = NA_integer_,
    tag = "",
    show_progress = TRUE
  )

  if (length(args) == 0L) return(out)

  for (a in args) {
    if (!grepl("^--", a)) stop("Bad arg: ", a)
    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
    key <- kv[1]
    val <- if (length(kv) > 1L) kv[2] else "TRUE"

    if (key == "n") out$n <- as.integer(val)
    else if (key == "times") out$times <- as.integer(val)
    else if (key == "base_seed") out$base_seed <- as.integer(val)
    else if (key == "nmulti") out$nmulti <- as.integer(val)
    else if (key == "lp_bases") out$lp_bases <- parse_csv(val)
    else if (key == "lp_degree") out$lp_degree <- val
    else if (key == "lp_bernstein.basis" || key == "lp_bernstein" || key == "lp_bernstein_basis") out$lp_bernstein.basis <- parse_csv_bool(val)
    else if (key == "nslaves" || key == "rslaves") out$nslaves <- as.integer(val)
    else if (key == "out_dir") out$out_dir <- val
    else if (key == "fi_tcp_iface") out$fi_tcp_iface <- val
    else if (key == "r_libs") out$r_libs <- val
    else if (key == "max_combos") out$max_combos <- as.integer(val)
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
  sprintf("%s_%s_%s_n%d_t%d_s%d_m%d_ns%d%s", backend, hash, ts, cfg$n, cfg$times, cfg$base_seed, cfg$nmulti, cfg$nslaves, tag)
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  cfg <- parse_args(args)

  argv <- commandArgs(trailingOnly = FALSE)
  script_arg <- argv[grep("^--file=", argv)][1]
  script_path <- if (!is.na(script_arg) && nzchar(script_arg)) {
    sub("^--file=", "", script_arg)
  } else {
    file.path(getwd(), "run_npreg_combos.R")
  }
  here <- dirname(normalizePath(script_path))

  bench_script <- file.path(here, "bench_npreg_param_nprmpi.R")
  if (!file.exists(bench_script)) stop("Missing benchmark script: ", bench_script)

  hash <- repo_hash(here)
  run_id <- mk_run_id(cfg, "npRmpi", hash)

  out_raw <- file.path(cfg$out_dir, paste0("npreg_combo_raw_", run_id, ".csv"))
  out_summary <- file.path(cfg$out_dir, paste0("npreg_combo_summary_", run_id, ".csv"))
  out_manifest <- file.path(cfg$out_dir, paste0("npreg_combo_manifest_", run_id, ".csv"))

  base_grid <- expand.grid(
    bwmethod = c("cv.ls", "cv.aic"),
    ckertype = c("gaussian", "epanechnikov"),
    np_tree = c(FALSE, TRUE),
    seed_policy = c("fixed", "varying"),
    stringsAsFactors = FALSE
  )

  combos_base <- cbind(
    regtype = rep(c("lc", "ll"), each = nrow(base_grid)),
    base_grid[rep(seq_len(nrow(base_grid)), times = 2L), , drop = FALSE],
    basis = NA_character_,
    degree = "",
    bernstein.basis = NA
  )

  combos_lp <- expand.grid(
    regtype = "lp",
    basis = cfg$lp_bases,
    degree = cfg$lp_degree,
    bernstein.basis = cfg$lp_bernstein.basis,
    bwmethod = c("cv.ls", "cv.aic"),
    ckertype = c("gaussian", "epanechnikov"),
    np_tree = c(FALSE, TRUE),
    seed_policy = c("fixed", "varying"),
    stringsAsFactors = FALSE
  )

  combos <- rbind(
    combos_base[, c("regtype", "basis", "degree", "bernstein.basis", "bwmethod", "ckertype", "np_tree", "seed_policy")],
    combos_lp[, c("regtype", "basis", "degree", "bernstein.basis", "bwmethod", "ckertype", "np_tree", "seed_policy")]
  )

  status <- vector("list", nrow(combos))
  combo_n <- if (is.na(cfg$max_combos)) nrow(combos) else min(nrow(combos), cfg$max_combos)

  for (i in seq_len(combo_n)) {
    row <- combos[i, ]

    if (cfg$show_progress) {
      cat(sprintf("[%02d/%02d] regtype=%s basis=%s degree=%s bern=%s bw=%s cker=%s tree=%s seed_policy=%s\n",
                  i, nrow(combos), row$regtype,
                  ifelse(is.na(row$basis), "", row$basis),
                  ifelse(nchar(row$degree) == 0L, "", row$degree),
                  ifelse(is.na(row$bernstein.basis), "", toupper(as.character(row$bernstein.basis))),
                  row$bwmethod, row$ckertype,
                  as.character(row$np_tree), row$seed_policy))
    }

    cmd_args <- c(
      bench_script,
      sprintf("--n=%d", cfg$n),
      sprintf("--times=%d", cfg$times),
      sprintf("--base_seed=%d", cfg$base_seed),
      sprintf("--nmulti=%d", cfg$nmulti),
      sprintf("--nslaves=%d", cfg$nslaves),
      sprintf("--regtype=%s", row$regtype),
      sprintf("--basis=%s", ifelse(is.na(row$basis), "glp", row$basis)),
      sprintf("--degree=%s", ifelse(nchar(row$degree) == 0L, "2,2", row$degree)),
      sprintf("--bernstein.basis=%s", ifelse(is.na(row$bernstein.basis), "FALSE", toupper(as.character(row$bernstein.basis)))),
      sprintf("--bwmethod=%s", row$bwmethod),
      sprintf("--ckertype=%s", row$ckertype),
      sprintf("--np_tree=%s", toupper(as.character(row$np_tree))),
      sprintf("--seed_policy=%s", row$seed_policy),
      sprintf("--out_raw=%s", out_raw),
      sprintf("--out_summary=%s", out_summary),
      "--show_progress=FALSE"
    )

    combo_log <- file.path(cfg$out_dir, sprintf("npreg_combo_%s_%02d.log", run_id, i))
    env_vec <- c(
      if (nzchar(cfg$r_libs)) sprintf("R_LIBS=%s", cfg$r_libs) else NULL,
      if (nzchar(cfg$fi_tcp_iface)) sprintf("FI_TCP_IFACE=%s", cfg$fi_tcp_iface) else NULL
    )
    rc <- system2("env", c(env_vec, "script", "-q", "/dev/null", "Rscript", cmd_args), stdout = combo_log, stderr = combo_log)

    status[[i]] <- data.frame(
      combo_id = i,
      regtype = row$regtype,
      basis = row$basis,
      degree = row$degree,
      bernstein.basis = row$bernstein.basis,
      bwmethod = row$bwmethod,
      ckertype = row$ckertype,
      np_tree = row$np_tree,
      seed_policy = row$seed_policy,
      log = combo_log,
      rc = rc,
      stringsAsFactors = FALSE
    )
  }

  manifest <- do.call(rbind, status)
  write.csv(manifest, out_manifest, row.names = FALSE)

  cat("raw_csv=", out_raw, "\n", sep = "")
  cat("summary_csv=", out_summary, "\n", sep = "")
  cat("manifest_csv=", out_manifest, "\n", sep = "")

  if (any(manifest$rc != 0L)) {
    stop("One or more combos failed. See manifest: ", out_manifest)
  }
}

if (sys.nframe() == 0) main()
