#!/usr/bin/env Rscript

parse_args <- function(args) {
  out <- list(
    n_values = c(250L, 500L),
    times = 2L,
    base_seed = 42L,
    out_dir = file.path("/tmp", paste0("np_hotspot_wp0_pre_", format(Sys.time(), "%Y%m%d_%H%M%S"))),
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
    else if (key == "base_seed") out$base_seed <- as.integer(val)
    else if (key == "out_dir") out$out_dir <- val
    else if (key == "show_progress") out$show_progress <- as.logical(val)
    else stop("Unknown arg: ", key)
  }

  if (length(out$n_values) == 0L || anyNA(out$n_values) || any(out$n_values < 50L))
    stop("n_values must be comma-separated integers >= 50")
  if (!is.finite(out$times) || out$times < 1L)
    stop("times must be >= 1")

  out
}

quiet_eval <- function(expr) {
  tf <- tempfile(pattern = "np_hotspot_quiet_")
  con <- file(tf, open = "wt")
  sink(con)
  on.exit({
    sink()
    close(con)
    unlink(tf)
  }, add = TRUE)
  force(expr)
}

first_num <- function(...) {
  xs <- list(...)
  for (x in xs) {
    if (!is.null(x) && length(x) > 0L) {
      v <- suppressWarnings(as.numeric(unlist(x)))
      v <- v[is.finite(v)]
      if (length(v) > 0L) return(v[[1L]])
    }
  }
  NA_real_
}

sig_head <- function(x, k = 20L) {
  v <- suppressWarnings(as.numeric(unlist(x)))
  v <- v[is.finite(v)]
  if (length(v) == 0L) return("")
  v <- v[seq_len(min(k, length(v)))]
  paste(signif(v, 10), collapse = ";")
}

bw_metrics <- function(bw) {
  if (is.null(bw)) {
    return(list(
      num_feval = NA_real_,
      fval = NA_real_,
      ifval = NA_real_,
      obj_hist_len = NA_integer_,
      obj_hist_sig = "",
      bw_sig = ""
    ))
  }

  hist_obj <- bw$eval.history
  hist_vals <- suppressWarnings(as.numeric(unlist(hist_obj)))
  hist_vals <- hist_vals[is.finite(hist_vals)]

  list(
    num_feval = first_num(bw$num.feval, bw$num.fval, bw$num.fevals),
    fval = first_num(bw$fval),
    ifval = first_num(bw$ifval),
    obj_hist_len = if (is.null(hist_obj)) NA_integer_ else as.integer(length(hist_vals)),
    obj_hist_sig = sig_head(hist_vals, k = 20L),
    bw_sig = sig_head(bw$bw, k = 20L)
  )
}

output_metrics <- function(out) {
  vals <- suppressWarnings(as.numeric(unlist(out)))
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0L) {
    return(list(out_len = 0L, out_mean = NA_real_, out_sd = NA_real_, out_sig = ""))
  }
  list(
    out_len = as.integer(length(vals)),
    out_mean = mean(vals),
    out_sd = stats::sd(vals),
    out_sig = sig_head(vals, k = 20L)
  )
}

mk_reg_data <- function(n, seed) {
  set.seed(seed)
  x1 <- runif(n)
  x2 <- runif(n)
  z1 <- factor(rbinom(n, 3, 0.5))
  z2 <- ordered(rbinom(n, 3, 0.5))
  y <- sin(2 * pi * x1) + x2 + as.numeric(z1) + rnorm(n, sd = 0.25)
  d <- data.frame(y = y, x1 = x1, x2 = x2, z1 = z1, z2 = z2)
  list(train = d, eval = d[seq_len(min(100L, n)), c("x1", "x2", "z1", "z2")])
}

mk_scoef_data <- function(n, seed) {
  set.seed(seed)
  x <- runif(n)
  z <- runif(n, -2, 2)
  y <- x * exp(z) * (1 + rnorm(n, sd = 0.2))
  d <- data.frame(y = y, x = x, z = z)
  idx <- seq_len(min(100L, n))
  list(
    train = d,
    ex = d[idx, "x", drop = FALSE],
    ez = d[idx, "z", drop = FALSE]
  )
}

mk_index_data <- function(n, seed) {
  set.seed(seed)
  x1 <- runif(n, -1, 1)
  x2 <- runif(n, -1, 1)
  y <- x1 - x2 + rnorm(n)
  d <- data.frame(y = y, x1 = x1, x2 = x2)
  list(train = d, eval = d[seq_len(min(100L, n)), c("x1", "x2")])
}

scenario_definitions <- function() {
  list(
    list(
      id = "npreg_bw_ll", function_name = "npreg", path = "bw", regtype = "ll", basis_degree = "-",
      runner = function(n, seed) {
        dat <- mk_reg_data(n, seed)
        bw <- npregbw(y ~ x1 + x2 + z1 + z2, data = dat$train, regtype = "ll", bwmethod = "cv.ls", nmulti = 1)
        list(bw = bw, out = bw$bw)
      }
    ),
    list(
      id = "npreg_bw_lp", function_name = "npreg", path = "bw", regtype = "lp", basis_degree = "glp/2,2",
      runner = function(n, seed) {
        dat <- mk_reg_data(n, seed)
        bw <- npregbw(
          y ~ x1 + x2 + z1 + z2, data = dat$train, regtype = "lp",
          basis = "glp", degree = c(2L, 2L), bwmethod = "cv.ls", nmulti = 1
        )
        list(bw = bw, out = bw$bw)
      }
    ),
    list(
      id = "npreghat_matrix_ll", function_name = "npreghat", path = "hat_matrix", regtype = "ll", basis_degree = "-",
      runner = function(n, seed) {
        dat <- mk_reg_data(n, seed)
        bw <- npregbw(y ~ x1 + x2 + z1 + z2, data = dat$train, regtype = "ll", bwmethod = "cv.ls", nmulti = 1)
        H <- npreghat(
          bws = bw,
          txdat = dat$train[, c("x1", "x2", "z1", "z2")],
          tydat = dat$train$y,
          exdat = dat$eval,
          output = "matrix"
        )
        list(bw = bw, out = H)
      }
    ),
    list(
      id = "npreghat_matrix_lp", function_name = "npreghat", path = "hat_matrix", regtype = "lp", basis_degree = "glp/2,2",
      runner = function(n, seed) {
        dat <- mk_reg_data(n, seed)
        bw <- npregbw(
          y ~ x1 + x2 + z1 + z2, data = dat$train, regtype = "lp",
          basis = "glp", degree = c(2L, 2L), bwmethod = "cv.ls", nmulti = 1
        )
        H <- npreghat(
          bws = bw,
          txdat = dat$train[, c("x1", "x2", "z1", "z2")],
          tydat = dat$train$y,
          exdat = dat$eval,
          output = "matrix"
        )
        list(bw = bw, out = H)
      }
    ),
    list(
      id = "npscoef_bw_ll", function_name = "npscoef", path = "bw", regtype = "ll", basis_degree = "-",
      runner = function(n, seed) {
        dat <- mk_scoef_data(n, seed)
        bw <- npscoefbw(y ~ x | z, data = dat$train, regtype = "ll", bwmethod = "cv.ls", nmulti = 1)
        list(bw = bw, out = bw$bw)
      }
    ),
    list(
      id = "npscoef_bw_lp", function_name = "npscoef", path = "bw", regtype = "lp", basis_degree = "glp/2",
      runner = function(n, seed) {
        dat <- mk_scoef_data(n, seed)
        bw <- npscoefbw(y ~ x | z, data = dat$train, regtype = "lp", basis = "glp", degree = 2L, bwmethod = "cv.ls", nmulti = 1)
        list(bw = bw, out = bw$bw)
      }
    ),
    list(
      id = "npscoefhat_apply_ll", function_name = "npscoefhat", path = "hat_apply", regtype = "ll", basis_degree = "-",
      runner = function(n, seed) {
        dat <- mk_scoef_data(n, seed)
        bw <- npscoefbw(y ~ x | z, data = dat$train, regtype = "ll", bwmethod = "cv.ls", nmulti = 1)
        out <- npscoefhat(
          bws = bw,
          txdat = dat$train["x"],
          tzdat = dat$train["z"],
          exdat = dat$ex,
          ezdat = dat$ez,
          output = "apply",
          y = dat$train$y
        )
        list(bw = bw, out = out)
      }
    ),
    list(
      id = "npscoefhat_apply_lp", function_name = "npscoefhat", path = "hat_apply", regtype = "lp", basis_degree = "glp/2",
      runner = function(n, seed) {
        dat <- mk_scoef_data(n, seed)
        bw <- npscoefbw(y ~ x | z, data = dat$train, regtype = "lp", basis = "glp", degree = 2L, bwmethod = "cv.ls", nmulti = 1)
        out <- npscoefhat(
          bws = bw,
          txdat = dat$train["x"],
          tzdat = dat$train["z"],
          exdat = dat$ex,
          ezdat = dat$ez,
          output = "apply",
          y = dat$train$y
        )
        list(bw = bw, out = out)
      }
    ),
    list(
      id = "npindex_bw_ll", function_name = "npindex", path = "bw", regtype = "ll", basis_degree = "-",
      runner = function(n, seed) {
        dat <- mk_index_data(n, seed)
        bw <- npindexbw(y ~ x1 + x2, data = dat$train, method = "ichimura", regtype = "ll", nmulti = 1)
        list(bw = bw, out = bw$bw)
      }
    ),
    list(
      id = "npindex_bw_lp", function_name = "npindex", path = "bw", regtype = "lp", basis_degree = "glp/2",
      runner = function(n, seed) {
        dat <- mk_index_data(n, seed)
        bw <- npindexbw(
          y ~ x1 + x2, data = dat$train, method = "ichimura", regtype = "lp",
          basis = "glp", degree = 2L, nmulti = 1
        )
        list(bw = bw, out = bw$bw)
      }
    ),
    list(
      id = "npindexhat_apply_ll", function_name = "npindexhat", path = "hat_apply", regtype = "ll", basis_degree = "-",
      runner = function(n, seed) {
        dat <- mk_index_data(n, seed)
        bw <- npindexbw(y ~ x1 + x2, data = dat$train, method = "ichimura", regtype = "ll", nmulti = 1)
        out <- npindexhat(
          bws = bw,
          txdat = dat$train[, c("x1", "x2")],
          exdat = dat$eval,
          output = "apply",
          y = dat$train$y
        )
        list(bw = bw, out = out)
      }
    ),
    list(
      id = "npindexhat_apply_lp", function_name = "npindexhat", path = "hat_apply", regtype = "lp", basis_degree = "glp/2",
      runner = function(n, seed) {
        dat <- mk_index_data(n, seed)
        bw <- npindexbw(
          y ~ x1 + x2, data = dat$train, method = "ichimura", regtype = "lp",
          basis = "glp", degree = 2L, nmulti = 1
        )
        out <- npindexhat(
          bws = bw,
          txdat = dat$train[, c("x1", "x2")],
          exdat = dat$eval,
          output = "apply",
          y = dat$train$y
        )
        list(bw = bw, out = out)
      }
    )
  )
}

run_one <- function(sc, n, seed_policy, iter, base_seed, sc_index) {
  seed <- if (seed_policy == "fixed") base_seed else base_seed + (sc_index * 1000L) + (n * 10L) + iter

  gc(reset = TRUE)
  t0 <- proc.time()[["elapsed"]]
  got <- tryCatch(
    quiet_eval(sc$runner(n = n, seed = seed)),
    error = function(e) e
  )
  elapsed <- proc.time()[["elapsed"]] - t0
  g1 <- gc()

  if (inherits(got, "error")) {
    return(data.frame(
      scenario_id = sc$id,
      function_name = sc$function_name,
      path = sc$path,
      regtype = sc$regtype,
      basis_degree = sc$basis_degree,
      n = n,
      seed_policy = seed_policy,
      iter = iter,
      seed = seed,
      elapsed_sec = elapsed,
      num_feval = NA_real_,
      fval = NA_real_,
      ifval = NA_real_,
      obj_hist_len = NA_integer_,
      obj_hist_sig = "",
      bw_sig = "",
      out_len = NA_integer_,
      out_mean = NA_real_,
      out_sd = NA_real_,
      out_sig = "",
      gc_ncells_peak_mb = as.numeric(g1["Ncells", 7L]),
      gc_vcells_peak_mb = as.numeric(g1["Vcells", 7L]),
      ok = FALSE,
      error = conditionMessage(got),
      stringsAsFactors = FALSE
    ))
  }

  bwm <- bw_metrics(got$bw)
  om <- output_metrics(got$out)

  data.frame(
    scenario_id = sc$id,
    function_name = sc$function_name,
    path = sc$path,
    regtype = sc$regtype,
    basis_degree = sc$basis_degree,
    n = n,
    seed_policy = seed_policy,
    iter = iter,
    seed = seed,
    elapsed_sec = elapsed,
    num_feval = bwm$num_feval,
    fval = bwm$fval,
    ifval = bwm$ifval,
    obj_hist_len = bwm$obj_hist_len,
    obj_hist_sig = bwm$obj_hist_sig,
    bw_sig = bwm$bw_sig,
    out_len = om$out_len,
    out_mean = om$out_mean,
    out_sd = om$out_sd,
    out_sig = om$out_sig,
    gc_ncells_peak_mb = as.numeric(g1["Ncells", 7L]),
    gc_vcells_peak_mb = as.numeric(g1["Vcells", 7L]),
    ok = TRUE,
    error = "",
    stringsAsFactors = FALSE
  )
}

summarize_group <- function(df) {
  ok <- df[df$ok, , drop = FALSE]
  if (nrow(ok) == 0L) {
    return(data.frame(
      scenario_id = df$scenario_id[1],
      function_name = df$function_name[1],
      path = df$path[1],
      regtype = df$regtype[1],
      basis_degree = df$basis_degree[1],
      n = df$n[1],
      seed_policy = df$seed_policy[1],
      runs = nrow(df),
      ok_runs = 0L,
      mean_elapsed_sec = NA_real_,
      median_elapsed_sec = NA_real_,
      num_feval_unique = NA_integer_,
      fval_min = NA_real_,
      fval_max = NA_real_,
      fval_range = NA_real_,
      obj_hist_len_unique = NA_integer_,
      bw_sig_unique = NA_integer_,
      mean_gc_vcells_peak_mb = NA_real_,
      status = "FAIL",
      stringsAsFactors = FALSE
    ))
  }

  fvals <- ok$fval[is.finite(ok$fval)]
  numf <- ok$num_feval[is.finite(ok$num_feval)]
  objlen <- ok$obj_hist_len[is.finite(ok$obj_hist_len)]
  bws <- ok$bw_sig[nzchar(ok$bw_sig)]

  data.frame(
    scenario_id = ok$scenario_id[1],
    function_name = ok$function_name[1],
    path = ok$path[1],
    regtype = ok$regtype[1],
    basis_degree = ok$basis_degree[1],
    n = ok$n[1],
    seed_policy = ok$seed_policy[1],
    runs = nrow(df),
    ok_runs = nrow(ok),
    mean_elapsed_sec = mean(ok$elapsed_sec),
    median_elapsed_sec = stats::median(ok$elapsed_sec),
    num_feval_unique = if (length(numf)) length(unique(numf)) else NA_integer_,
    fval_min = if (length(fvals)) min(fvals) else NA_real_,
    fval_max = if (length(fvals)) max(fvals) else NA_real_,
    fval_range = if (length(fvals)) max(fvals) - min(fvals) else NA_real_,
    obj_hist_len_unique = if (length(objlen)) length(unique(objlen)) else NA_integer_,
    bw_sig_unique = if (length(bws)) length(unique(bws)) else NA_integer_,
    mean_gc_vcells_peak_mb = mean(ok$gc_vcells_peak_mb),
    status = if (nrow(ok) == nrow(df)) "OK" else "PARTIAL",
    stringsAsFactors = FALSE
  )
}

build_prepost_scaffold <- function(sumdf) {
  data.frame(
    `function` = sumdf$function_name,
    path = sumdf$path,
    regtype = sumdf$regtype,
    basis_degree = sumdf$basis_degree,
    n = sumdf$n,
    seed_policy = sumdf$seed_policy,
    num_feval_pre = NA_real_,
    num_feval_post = NA_real_,
    fval_pre = NA_real_,
    fval_post = NA_real_,
    parity_status = "BASELINE_CAPTURED",
    pre_mean_sec = sumdf$mean_elapsed_sec,
    post_mean_sec = NA_real_,
    ratio_post_pre = NA_real_,
    delta_pct = NA_real_,
    stringsAsFactors = FALSE
  )
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  cfg <- parse_args(args)
  suppressPackageStartupMessages(library(np))
  options(np.messages = FALSE, np.tree = FALSE)

  dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
  scenarios <- scenario_definitions()

  if (cfg$show_progress) {
    cat("Running hotspot WP0/WP1 pre-baseline in np-master\n")
    cat("out_dir:", cfg$out_dir, "\n")
    cat("n_values:", paste(cfg$n_values, collapse = ","), " times:", cfg$times, "\n")
  }

  rows <- list()
  idx <- 0L
  total_jobs <- length(scenarios) * length(cfg$n_values) * 2L * cfg$times
  job <- 0L

  for (s in seq_along(scenarios)) {
    sc <- scenarios[[s]]
    for (n in cfg$n_values) {
      for (seed_policy in c("fixed", "varying")) {
        for (iter in seq_len(cfg$times)) {
          job <- job + 1L
          if (cfg$show_progress) {
            cat(sprintf("[%d/%d] %s n=%d seed_policy=%s iter=%d\n", job, total_jobs, sc$id, n, seed_policy, iter))
          }
          idx <- idx + 1L
          rows[[idx]] <- run_one(
            sc = sc,
            n = as.integer(n),
            seed_policy = seed_policy,
            iter = as.integer(iter),
            base_seed = cfg$base_seed,
            sc_index = as.integer(s)
          )
        }
      }
    }
  }

  raw <- do.call(rbind, rows)
  raw_path <- file.path(cfg$out_dir, "wp0_wp1_pre_raw.csv")
  write.csv(raw, raw_path, row.names = FALSE)

  split_key <- interaction(raw$scenario_id, raw$n, raw$seed_policy, drop = TRUE)
  sums <- do.call(rbind, lapply(split(raw, split_key), summarize_group))
  sums <- sums[order(sums$function_name, sums$path, sums$regtype, sums$n, sums$seed_policy), ]
  sum_path <- file.path(cfg$out_dir, "wp0_wp1_pre_summary.csv")
  write.csv(sums, sum_path, row.names = FALSE)

  scaffold <- build_prepost_scaffold(sums)
  scaffold_path <- file.path(cfg$out_dir, "wp0_wp1_prepost_scaffold.csv")
  write.csv(scaffold, scaffold_path, row.names = FALSE)

  meta <- data.frame(
    generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    n_values = paste(cfg$n_values, collapse = ","),
    times = cfg$times,
    base_seed = cfg$base_seed,
    stringsAsFactors = FALSE
  )
  write.csv(meta, file.path(cfg$out_dir, "wp0_wp1_pre_meta.csv"), row.names = FALSE)

  if (cfg$show_progress) {
    cat("\nArtifacts:\n")
    cat("  raw:", raw_path, "\n")
    cat("  summary:", sum_path, "\n")
    cat("  scaffold:", scaffold_path, "\n")
    cat("\nSummary preview:\n")
    print(utils::head(sums, 12))
  }
}

if (sys.nframe() == 0L) main()
