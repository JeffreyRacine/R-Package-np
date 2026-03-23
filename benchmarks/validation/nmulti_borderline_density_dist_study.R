#!/usr/bin/env Rscript

parse_args <- function(args) {
  out <- list(
    n_grid = c(125L, 250L),
    seeds = 101:106,
    nmulti_grid = c(1L, 2L, 5L),
    out_dir = sprintf("/tmp/nmulti_borderline_density_dist_%s", format(Sys.time(), "%Y%m%d_%H%M%S")),
    show_progress = TRUE
  )

  parse_int_csv <- function(x, name) {
    vals <- as.integer(strsplit(x, ",", fixed = TRUE)[[1L]])
    if (length(vals) == 0L || anyNA(vals)) stop(name, " must be a comma-separated integer list")
    vals
  }

  if (length(args) == 0L) return(out)

  for (a in args) {
    if (!grepl("^--", a)) stop("Bad arg: ", a)
    raw <- sub("^--", "", a)
    key <- sub("=.*$", "", raw)
    val <- if (grepl("=", raw, fixed = TRUE)) sub("^[^=]*=", "", raw) else "TRUE"

    if (key == "n_grid") out$n_grid <- parse_int_csv(val, "n_grid")
    else if (key == "seeds") out$seeds <- parse_int_csv(val, "seeds")
    else if (key == "nmulti_grid") out$nmulti_grid <- parse_int_csv(val, "nmulti_grid")
    else if (key == "out_dir") out$out_dir <- val
    else if (key == "show_progress") out$show_progress <- as.logical(val)
    else stop("Unknown arg: ", key)
  }

  out$n_grid <- sort(unique(out$n_grid))
  out$seeds <- sort(unique(out$seeds))
  out$nmulti_grid <- sort(unique(out$nmulti_grid))
  if (!5L %in% out$nmulti_grid) stop("nmulti_grid must include 5 as the reference")
  if (any(out$nmulti_grid < 1L)) stop("nmulti_grid must be >= 1")
  out
}

load_backend <- function() {
  suppressPackageStartupMessages(library(np))
  options(np.messages = FALSE)
}

scenario_catalog <- function() {
  list(
    list(id = "npcdens_ll_dgp1", family = "npcdens", regtype = "ll", degree = 1L, dgp = "dgp1"),
    list(id = "npcdens_lp2_dgp1", family = "npcdens", regtype = "lp", degree = 2L, dgp = "dgp1"),
    list(id = "npcdist_ll_dgp1", family = "npcdist", regtype = "ll", degree = 1L, dgp = "dgp1"),
    list(id = "npcdist_lp2_dgp1", family = "npcdist", regtype = "lp", degree = 2L, dgp = "dgp1"),
    list(id = "npcdens_ll_dgp2", family = "npcdens", regtype = "ll", degree = 1L, dgp = "dgp2"),
    list(id = "npcdens_lp2_dgp2", family = "npcdens", regtype = "lp", degree = 2L, dgp = "dgp2"),
    list(id = "npcdist_ll_dgp2", family = "npcdist", regtype = "ll", degree = 1L, dgp = "dgp2"),
    list(id = "npcdist_lp2_dgp2", family = "npcdist", regtype = "lp", degree = 2L, dgp = "dgp2")
  )
}

sample_z <- function(n) {
  factor(sample(c("a", "b", "c"), n, replace = TRUE, prob = c(0.35, 0.40, 0.25)),
         levels = c("a", "b", "c"))
}

truth_dgp1_mu <- function(x, z) {
  z_eff <- c(a = -0.30, b = 0.20, c = 0.70)[as.character(z)]
  0.8 * sin(2 * pi * x) + 0.35 * cos(pi * x) + z_eff
}

truth_dgp1_sigma <- function(x, z) {
  z_scale <- c(a = 1.0, b = 1.15, c = 0.90)[as.character(z)]
  (0.18 + 0.12 * x) * z_scale
}

truth_dgp2_params <- function(x, z) {
  z_char <- as.character(z)
  w <- ifelse(z_char == "a", 0.75, ifelse(z_char == "b", 0.55, 0.35))
  mu1 <- -0.45 + 0.9 * sin(2 * pi * x) + ifelse(z_char == "c", 0.25, 0.0)
  mu2 <- 0.75 + 0.6 * cos(pi * x) + ifelse(z_char == "a", 0.15, -0.10)
  sd1 <- 0.16 + 0.10 * x
  sd2 <- 0.26 + 0.08 * (1 - x)
  list(w = w, mu1 = mu1, mu2 = mu2, sd1 = sd1, sd2 = sd2)
}

generate_data <- function(seed, n, dgp) {
  set.seed(seed)
  x <- runif(n)
  z <- sample_z(n)

  if (identical(dgp, "dgp1")) {
    mu <- truth_dgp1_mu(x, z)
    sigma <- truth_dgp1_sigma(x, z)
    y <- mu + sigma * rnorm(n)
  } else {
    pars <- truth_dgp2_params(x, z)
    comp <- rbinom(n, size = 1L, prob = pars$w)
    y <- ifelse(
      comp == 1L,
      pars$mu1 + pars$sd1 * rnorm(n),
      pars$mu2 + pars$sd2 * rnorm(n)
    )
  }

  data.frame(y = y, x = x, z = z)
}

build_eval_data <- function() {
  anchors <- data.frame(
    x = c(0.10, 0.25, 0.45, 0.60, 0.80, 0.90),
    z = factor(c("a", "b", "c", "a", "b", "c"), levels = c("a", "b", "c")),
    slice_id = seq_len(6L)
  )
  y_grid <- seq(-3.5, 4.5, length.out = 101L)
  rows <- vector("list", nrow(anchors))
  offset <- 0L
  for (i in seq_len(nrow(anchors))) {
    cur <- anchors[rep(i, length(y_grid)), , drop = FALSE]
    cur$y <- y_grid
    cur$eval_id <- seq_len(length(y_grid)) + offset
    offset <- offset + length(y_grid)
    rows[[i]] <- cur
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

truth_density <- function(eval_df, dgp) {
  x <- eval_df$x
  z <- eval_df$z
  y <- eval_df$y

  if (identical(dgp, "dgp1")) {
    mu <- truth_dgp1_mu(x, z)
    sigma <- truth_dgp1_sigma(x, z)
    return(dnorm(y, mean = mu, sd = sigma))
  }

  pars <- truth_dgp2_params(x, z)
  pars$w * dnorm(y, mean = pars$mu1, sd = pars$sd1) +
    (1 - pars$w) * dnorm(y, mean = pars$mu2, sd = pars$sd2)
}

truth_cdf <- function(eval_df, dgp) {
  x <- eval_df$x
  z <- eval_df$z
  y <- eval_df$y

  if (identical(dgp, "dgp1")) {
    mu <- truth_dgp1_mu(x, z)
    sigma <- truth_dgp1_sigma(x, z)
    return(pnorm(y, mean = mu, sd = sigma))
  }

  pars <- truth_dgp2_params(x, z)
  pars$w * pnorm(y, mean = pars$mu1, sd = pars$sd1) +
    (1 - pars$w) * pnorm(y, mean = pars$mu2, sd = pars$sd2)
}

trapz <- function(x, y) {
  n <- length(x)
  if (n < 2L) return(NA_real_)
  sum((x[-1L] - x[-n]) * (y[-1L] + y[-n]) / 2)
}

compute_surface_metric <- function(eval_df, pred, truth) {
  err <- pred - truth
  slices <- split(seq_len(nrow(eval_df)), eval_df$slice_id)
  iae <- numeric(length(slices))
  ise <- numeric(length(slices))
  for (i in seq_along(slices)) {
    idx <- slices[[i]]
    y <- eval_df$y[idx]
    iae[i] <- trapz(y, abs(err[idx]))
    ise[i] <- trapz(y, err[idx]^2)
  }
  list(
    primary = mean(iae),
    iae = mean(iae),
    ise = mean(ise),
    rmse = sqrt(mean(err^2))
  )
}

encode_numeric <- function(x) paste(signif(as.numeric(x), 10L), collapse = ";")

extract_bw_values <- function(bw) {
  if (!is.null(bw$bandwidth) && length(bw$bandwidth) > 0L) return(as.numeric(unlist(bw$bandwidth)))
  c(as.numeric(bw$ybw), as.numeric(bw$xbw))
}

count_history <- function(x) {
  if (is.null(x)) return(NA_integer_)
  nr <- nrow(x)
  if (!is.null(nr) && length(nr) > 0L) return(as.integer(nr)[1L])
  as.integer(length(x))[1L]
}

first_num_or_na <- function(x) {
  if (is.null(x) || length(x) == 0L) return(NA_real_)
  as.numeric(x)[1L]
}

fit_one <- function(spec, dat, nmulti, eval_df) {
  bw_fun <- if (identical(spec$family, "npcdens")) npcdensbw else npcdistbw
  fit_fun <- if (identical(spec$family, "npcdens")) npcdens else npcdist
  bw_args <- list(
    formula = y ~ x + z,
    data = dat,
    regtype = spec$regtype,
    bwtype = "fixed",
    bwmethod = if (identical(spec$family, "npcdens")) "cv.ml" else "cv.ls",
    nmulti = nmulti
  )
  if (identical(spec$regtype, "lp")) {
    bw_args$degree <- spec$degree
    bw_args$basis <- "glp"
    bw_args$bernstein.basis <- FALSE
  }

  t0 <- proc.time()[["elapsed"]]
  bw <- do.call(bw_fun, bw_args)
  elapsed_bw <- proc.time()[["elapsed"]] - t0
  fit <- fit_fun(bws = bw, data = dat)

  pred <- as.numeric(predict(fit, newdata = eval_df[, c("y", "x", "z"), drop = FALSE]))
  truth <- if (identical(spec$family, "npcdens")) truth_density(eval_df, spec$dgp) else truth_cdf(eval_df, spec$dgp)
  surf <- compute_surface_metric(eval_df, pred, truth)

  eval_out <- data.frame(
    scenario_id = spec$id,
    family = spec$family,
    dgp = spec$dgp,
    regtype = spec$regtype,
    degree = spec$degree,
    n = nrow(dat),
    seed = NA_integer_,
    nmulti = nmulti,
    eval_id = eval_df$eval_id,
    slice_id = eval_df$slice_id,
    pred = pred,
    truth = truth,
    stringsAsFactors = FALSE
  )

  list(
    bw = bw,
    elapsed_bw = elapsed_bw,
    surface = surf,
    eval = eval_out
  )
}

run_one <- function(spec, n, seed, nmulti, eval_df) {
  dat <- generate_data(seed = seed, n = n, dgp = spec$dgp)
  res <- fit_one(spec = spec, dat = dat, nmulti = nmulti, eval_df = eval_df)
  res$eval$seed <- seed
  list(
    run = data.frame(
      scenario_id = spec$id,
      family = spec$family,
      dgp = spec$dgp,
      regtype = spec$regtype,
      degree = spec$degree,
      n = n,
      seed = seed,
      nmulti = nmulti,
      elapsed_bw = res$elapsed_bw,
      fval = first_num_or_na(res$bw$fval),
      ifval = first_num_or_na(res$bw$ifval),
      num_fval = first_num_or_na(if (!is.null(res$bw$num.fval)) res$bw$num.fval else res$bw$num.feval),
      eval_count = count_history(res$bw$eval.history),
      invalid_count = count_history(res$bw$invalid.history),
      bandwidth = encode_numeric(extract_bw_values(res$bw)),
      surface_primary = res$surface$primary,
      surface_iae = res$surface$iae,
      surface_ise = res$surface$ise,
      surface_rmse = res$surface$rmse,
      stringsAsFactors = FALSE
    ),
    eval = res$eval
  )
}

bandwidth_l2_diff <- function(a, b) {
  va <- as.numeric(strsplit(a, ";", fixed = TRUE)[[1L]])
  vb <- as.numeric(strsplit(b, ";", fixed = TRUE)[[1L]])
  sqrt(sum((va - vb)^2))
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  median(x)
}

safe_quantile <- function(x, prob) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  as.numeric(stats::quantile(x, probs = prob, names = FALSE))
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  cfg <- parse_args(args)
  dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
  load_backend()

  scenarios <- scenario_catalog()
  eval_df <- build_eval_data()
  total_jobs <- length(scenarios) * length(cfg$n_grid) * length(cfg$seeds) * length(cfg$nmulti_grid)
  job_idx <- 0L

  raw_runs <- list()
  raw_eval <- list()
  run_idx <- 0L
  eval_idx <- 0L

  if (cfg$show_progress) {
    cat("Running borderline density/distribution study\n")
    cat("out_dir=", cfg$out_dir, "\n", sep = "")
  }

  for (spec in scenarios) {
    for (n in cfg$n_grid) {
      for (seed in cfg$seeds) {
        for (nmulti in cfg$nmulti_grid) {
          job_idx <- job_idx + 1L
          if (cfg$show_progress) cat(sprintf("[%d/%d] %s n=%d seed=%d nmulti=%d\n", job_idx, total_jobs, spec$id, n, seed, nmulti))
          res <- run_one(spec = spec, n = n, seed = seed, nmulti = nmulti, eval_df = eval_df)
          run_idx <- run_idx + 1L
          raw_runs[[run_idx]] <- res$run
          eval_idx <- eval_idx + 1L
          raw_eval[[eval_idx]] <- res$eval
        }
      }
    }
  }

  raw_runs_df <- do.call(rbind, raw_runs)
  raw_eval_df <- do.call(rbind, raw_eval)

  key <- interaction(raw_runs_df$scenario_id, raw_runs_df$n, raw_runs_df$seed, drop = TRUE)
  best_primary <- tapply(raw_runs_df$surface_primary, key, min, na.rm = TRUE)
  raw_runs_df$best_primary <- as.numeric(best_primary[key])
  raw_runs_df$within_1pct_best <- raw_runs_df$surface_primary <= 1.01 * raw_runs_df$best_primary

  ref_runs <- raw_runs_df[raw_runs_df$nmulti == 5L, , drop = FALSE]
  names(ref_runs)[names(ref_runs) %in% c("elapsed_bw", "fval", "bandwidth", "surface_primary")] <-
    c("elapsed_bw_ref5", "fval_ref5", "bandwidth_ref5", "surface_primary_ref5")

  merged <- merge(
    raw_runs_df,
    ref_runs[, c("scenario_id", "n", "seed", "elapsed_bw_ref5", "fval_ref5", "bandwidth_ref5", "surface_primary_ref5")],
    by = c("scenario_id", "n", "seed"),
    all.x = TRUE,
    sort = FALSE
  )
  merged$bw_speedup_vs_ref5 <- merged$elapsed_bw_ref5 / merged$elapsed_bw
  merged$fval_gap_vs_ref5 <- abs(merged$fval - merged$fval_ref5)
  merged$bandwidth_l2_vs_ref5 <- mapply(bandwidth_l2_diff, merged$bandwidth, merged$bandwidth_ref5)
  merged$surface_primary_delta_vs_ref5 <- merged$surface_primary - merged$surface_primary_ref5

  ref_eval <- raw_eval_df[raw_eval_df$nmulti == 5L, , drop = FALSE]
  names(ref_eval)[names(ref_eval) == "pred"] <- "pred_ref5"
  eval_merged <- merge(
    raw_eval_df,
    ref_eval[, c("scenario_id", "n", "seed", "eval_id", "pred_ref5")],
    by = c("scenario_id", "n", "seed", "eval_id"),
    all.x = TRUE,
    sort = FALSE
  )

  drift_groups <- split(eval_merged, interaction(eval_merged$scenario_id, eval_merged$n, eval_merged$seed, eval_merged$nmulti, drop = TRUE))
  drift_rows <- lapply(drift_groups, function(df) {
    data.frame(
      scenario_id = df$scenario_id[1L],
      family = df$family[1L],
      dgp = df$dgp[1L],
      regtype = df$regtype[1L],
      degree = df$degree[1L],
      n = df$n[1L],
      seed = df$seed[1L],
      nmulti = df$nmulti[1L],
      surface_drift_rmse_vs_ref5 = sqrt(mean((df$pred - df$pred_ref5)^2)),
      surface_drift_max_abs_vs_ref5 = max(abs(df$pred - df$pred_ref5)),
      stringsAsFactors = FALSE
    )
  })
  drift_df <- do.call(rbind, drift_rows)

  merged <- merge(
    merged,
    drift_df[, c("scenario_id", "n", "seed", "nmulti", "surface_drift_rmse_vs_ref5", "surface_drift_max_abs_vs_ref5")],
    by = c("scenario_id", "n", "seed", "nmulti"),
    all.x = TRUE,
    sort = FALSE
  )

  sum_groups <- split(merged, interaction(merged$scenario_id, merged$n, merged$nmulti, drop = TRUE))
  summary_rows <- lapply(sum_groups, function(df) {
    data.frame(
      scenario_id = df$scenario_id[1L],
      family = df$family[1L],
      dgp = df$dgp[1L],
      regtype = df$regtype[1L],
      degree = df$degree[1L],
      n = df$n[1L],
      nmulti = df$nmulti[1L],
      runs = nrow(df),
      median_elapsed_bw = safe_median(df$elapsed_bw),
      median_bw_speedup_vs_ref5 = safe_median(df$bw_speedup_vs_ref5),
      median_surface_primary = safe_median(df$surface_primary),
      median_surface_primary_delta_vs_ref5 = safe_median(df$surface_primary_delta_vs_ref5),
      median_surface_drift_rmse_vs_ref5 = safe_median(df$surface_drift_rmse_vs_ref5),
      p90_surface_drift_rmse_vs_ref5 = safe_quantile(df$surface_drift_rmse_vs_ref5, 0.90),
      median_fval_gap_vs_ref5 = safe_median(df$fval_gap_vs_ref5),
      median_bandwidth_l2_vs_ref5 = safe_median(df$bandwidth_l2_vs_ref5),
      share_within_1pct_best = mean(df$within_1pct_best),
      stringsAsFactors = FALSE
    )
  })
  summary_df <- do.call(rbind, summary_rows)

  runs_file <- file.path(cfg$out_dir, "nmulti_borderline_density_dist_runs.csv")
  eval_file <- file.path(cfg$out_dir, "nmulti_borderline_density_dist_eval.csv")
  drift_file <- file.path(cfg$out_dir, "nmulti_borderline_density_dist_drift.csv")
  summary_file <- file.path(cfg$out_dir, "nmulti_borderline_density_dist_summary.csv")
  report_file <- file.path(cfg$out_dir, "nmulti_borderline_density_dist_report.md")

  write.csv(raw_runs_df, runs_file, row.names = FALSE)
  write.csv(raw_eval_df, eval_file, row.names = FALSE)
  write.csv(drift_df, drift_file, row.names = FALSE)
  write.csv(summary_df, summary_file, row.names = FALSE)

  lines <- c(
    "# Borderline Density/Distribution nmulti Study",
    "",
    "Borderline routes only: `npcdens` / `npcdist`, `ll` and `lp(degree=2)`.",
    "Two mixed-data DGPs are used with one continuous and one unordered categorical regressor.",
    "`nmulti=5` is the practical reference.",
    ""
  )

  fams <- split(summary_df, summary_df$family)
  for (fam in names(fams)) {
    lines <- c(lines, sprintf("## %s", fam), "")
    by_scenario <- split(fams[[fam]], fams[[fam]]$scenario_id)
    for (sid in names(by_scenario)) {
      df <- by_scenario[[sid]]
      df <- df[order(df$n, df$nmulti), , drop = FALSE]
      lines <- c(lines, sprintf("### %s", sid), "")
      lines <- c(lines, "| n | nmulti | med bw sec | med bw speedup | med drift RMSE | p90 drift RMSE | share within 1% best |")
      lines <- c(lines, "|---:|---:|---:|---:|---:|---:|---:|")
      for (i in seq_len(nrow(df))) {
        lines <- c(lines, sprintf(
          "| %d | %d | %.3f | %.2f | %.6f | %.6f | %.2f |",
          df$n[i],
          df$nmulti[i],
          df$median_elapsed_bw[i],
          df$median_bw_speedup_vs_ref5[i],
          df$median_surface_drift_rmse_vs_ref5[i],
          df$p90_surface_drift_rmse_vs_ref5[i],
          df$share_within_1pct_best[i]
        ))
      }
      lines <- c(lines, "")
    }
  }
  writeLines(lines, con = report_file)

  if (cfg$show_progress) {
    cat("Runs:", runs_file, "\n")
    cat("Eval:", eval_file, "\n")
    cat("Drift:", drift_file, "\n")
    cat("Summary:", summary_file, "\n")
    cat("Report:", report_file, "\n")
    print(summary_df)
  }
}

if (sys.nframe() == 0L) main()
