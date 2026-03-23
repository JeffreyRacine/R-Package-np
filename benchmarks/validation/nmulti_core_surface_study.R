#!/usr/bin/env Rscript

parse_args <- function(args) {
  out <- list(
    n_grid = c(125L, 250L),
    seeds = 101:106,
    nmulti_grid = c(1L, 2L, 5L),
    out_dir = sprintf("/tmp/nmulti_core_surface_study_%s", format(Sys.time(), "%Y%m%d_%H%M%S")),
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

  out$n_grid <- sort(unique(as.integer(out$n_grid)))
  out$seeds <- sort(unique(as.integer(out$seeds)))
  out$nmulti_grid <- sort(unique(as.integer(out$nmulti_grid)))
  if (length(out$nmulti_grid) < 2L) stop("nmulti_grid must contain at least two values")
  if (any(out$nmulti_grid < 1L)) stop("nmulti_grid values must be >= 1")
  out
}

load_backend <- function() {
  suppressPackageStartupMessages(library(np))
  options(np.messages = FALSE)
}

scenario_catalog <- function() {
  list(
    list(id = "npreg_lc", family = "npreg", regtype = "lc", degree = NA_integer_),
    list(id = "npreg_ll", family = "npreg", regtype = "ll", degree = 1L),
    list(id = "npreg_lp_deg1", family = "npreg", regtype = "lp", degree = 1L),
    list(id = "npreg_lp_deg2", family = "npreg", regtype = "lp", degree = 2L),
    list(id = "npcdens_lc", family = "npcdens", regtype = "lc", degree = NA_integer_),
    list(id = "npcdens_ll", family = "npcdens", regtype = "ll", degree = 1L),
    list(id = "npcdens_lp_deg1", family = "npcdens", regtype = "lp", degree = 1L),
    list(id = "npcdens_lp_deg2", family = "npcdens", regtype = "lp", degree = 2L)
  )
}

generate_data <- function(seed, n) {
  set.seed(seed)
  x <- runif(n)
  z <- factor(sample(c("a", "b", "c"), n, replace = TRUE, prob = c(0.35, 0.40, 0.25)),
              levels = c("a", "b", "c"))
  mu <- truth_mu(data.frame(x = x, z = z))
  sigma <- truth_sigma(data.frame(x = x, z = z))
  y <- mu + sigma * rnorm(n)
  data.frame(y = y, x = x, z = z)
}

truth_mu <- function(df) {
  z_eff <- c(a = -0.25, b = 0.35, c = 0.80)[as.character(df$z)]
  0.80 * sin(2 * pi * df$x) + 0.40 * cos(pi * df$x) + z_eff
}

truth_sigma <- function(df) {
  0.20 + 0.15 * df$x
}

make_formula <- function() {
  y ~ x + z
}

build_regression_eval <- function() {
  data.frame(
    x = c(0.10, 0.30, 0.50, 0.70, 0.90, 0.20, 0.80, 0.55, 0.45),
    z = factor(c("a", "a", "a", "b", "b", "c", "c", "b", "c"), levels = c("a", "b", "c")),
    eval_id = seq_len(9L)
  )
}

build_density_eval <- function() {
  anchors <- data.frame(
    x = c(0.15, 0.45, 0.75, 0.25, 0.65, 0.85),
    z = factor(c("a", "a", "b", "c", "b", "c"), levels = c("a", "b", "c")),
    slice_id = seq_len(6L)
  )
  y_grid <- seq(-3.0, 4.5, length.out = 81L)
  rows <- vector("list", nrow(anchors))
  offset <- 0L
  for (i in seq_len(nrow(anchors))) {
    row_i <- anchors[rep(i, length(y_grid)), , drop = FALSE]
    row_i$y <- y_grid
    row_i$eval_id <- seq_len(length(y_grid)) + offset
    offset <- offset + length(y_grid)
    rows[[i]] <- row_i
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

trapz <- function(x, y) {
  n <- length(x)
  if (n < 2L) return(NA_real_)
  sum((x[-1L] - x[-n]) * (y[-1L] + y[-n]) / 2)
}

eval_truth <- function(spec, eval_df) {
  if (identical(spec$family, "npreg")) {
    return(truth_mu(eval_df))
  }
  mu <- truth_mu(eval_df)
  sigma <- truth_sigma(eval_df)
  dnorm(eval_df$y, mean = mu, sd = sigma)
}

compute_surface_metrics <- function(spec, eval_df, pred, truth) {
  err <- pred - truth
  if (identical(spec$family, "npreg")) {
    return(list(
      primary = sqrt(mean(err^2)),
      rmse = sqrt(mean(err^2)),
      mae = mean(abs(err)),
      iae = NA_real_,
      ise = NA_real_
    ))
  }

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
    rmse = sqrt(mean(err^2)),
    mae = mean(abs(err)),
    iae = mean(iae),
    ise = mean(ise)
  )
}

extract_bw_values <- function(bw) {
  vals <- NULL
  if (!is.null(bw$bandwidth) && length(bw$bandwidth) > 0L) {
    vals <- suppressWarnings(as.numeric(unlist(bw$bandwidth)))
  } else if (!is.null(bw$bw) && length(bw$bw) > 0L) {
    vals <- suppressWarnings(as.numeric(bw$bw))
  } else {
    vals <- c(
      suppressWarnings(as.numeric(bw$ybw)),
      suppressWarnings(as.numeric(bw$xbw))
    )
  }
  vals[is.finite(vals)]
}

encode_numeric <- function(x) {
  if (length(x) == 0L) return("")
  paste(signif(as.numeric(x), 10L), collapse = ";")
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
  bw_fun <- if (identical(spec$family, "npreg")) npregbw else npcdensbw
  fit_fun <- if (identical(spec$family, "npreg")) npreg else npcdens

  bw_args <- list(
    formula = make_formula(),
    data = dat,
    nmulti = nmulti,
    bwtype = "fixed",
    bwmethod = if (identical(spec$family, "npreg")) "cv.ls" else "cv.ml",
    regtype = spec$regtype
  )

  if (identical(spec$regtype, "lp")) {
    bw_args$basis <- "glp"
    bw_args$degree <- spec$degree
    bw_args$bernstein.basis <- FALSE
  }

  t0 <- proc.time()[["elapsed"]]
  bw <- do.call(bw_fun, bw_args)
  elapsed_bw <- proc.time()[["elapsed"]] - t0

  t1 <- proc.time()[["elapsed"]]
  fit <- fit_fun(bws = bw, data = dat)
  elapsed_fit <- proc.time()[["elapsed"]] - t1

  pred_data <- eval_df[, intersect(c("y", "x", "z"), names(eval_df)), drop = FALSE]
  pred <- as.numeric(predict(fit, newdata = pred_data))
  truth <- eval_truth(spec, eval_df)
  surf <- compute_surface_metrics(spec, eval_df, pred, truth)

  eval_out <- data.frame(
    scenario_id = spec$id,
    family = spec$family,
    regtype = spec$regtype,
    degree = ifelse(is.na(spec$degree), NA_integer_, spec$degree),
    n = nrow(dat),
    seed = NA_integer_,
    nmulti = nmulti,
    eval_id = eval_df$eval_id,
    slice_id = if ("slice_id" %in% names(eval_df)) eval_df$slice_id else NA_integer_,
    pred = pred,
    truth = truth,
    stringsAsFactors = FALSE
  )

  list(
    bw = bw,
    elapsed_bw = elapsed_bw,
    elapsed_fit = elapsed_fit,
    metrics = surf,
    eval = eval_out
  )
}

run_one <- function(spec, n, seed, nmulti, eval_df) {
  dat <- generate_data(seed = seed, n = n)
  res <- tryCatch(
    {
      fit <- fit_one(spec, dat, nmulti, eval_df)
      fit$eval$seed <- seed
      list(
        run = data.frame(
          scenario_id = spec$id,
          family = spec$family,
          regtype = spec$regtype,
          degree = ifelse(is.na(spec$degree), NA_integer_, spec$degree),
          n = n,
          seed = seed,
          nmulti = nmulti,
          elapsed_bw = fit$elapsed_bw,
          elapsed_fit = fit$elapsed_fit,
          elapsed_total = fit$elapsed_bw + fit$elapsed_fit,
          fval = first_num_or_na(fit$bw$fval),
          ifval = first_num_or_na(fit$bw$ifval),
          eval_count = count_history(fit$bw$eval.history),
          invalid_count = count_history(fit$bw$invalid.history),
          num_fval = first_num_or_na(if (!is.null(fit$bw$num.fval)) fit$bw$num.fval else fit$bw$num.feval),
          bandwidth = encode_numeric(extract_bw_values(fit$bw)),
          surface_primary = fit$metrics$primary,
          surface_rmse = fit$metrics$rmse,
          surface_mae = fit$metrics$mae,
          surface_iae = fit$metrics$iae,
          surface_ise = fit$metrics$ise,
          ok = TRUE,
          error = "",
          stringsAsFactors = FALSE
        ),
        eval = fit$eval
      )
    },
    error = function(e) {
      list(
        run = data.frame(
          scenario_id = spec$id,
          family = spec$family,
          regtype = spec$regtype,
          degree = ifelse(is.na(spec$degree), NA_integer_, spec$degree),
          n = n,
          seed = seed,
          nmulti = nmulti,
          elapsed_bw = NA_real_,
          elapsed_fit = NA_real_,
          elapsed_total = NA_real_,
          fval = NA_real_,
          ifval = NA_real_,
          eval_count = NA_integer_,
          invalid_count = NA_integer_,
          num_fval = NA_real_,
          bandwidth = "",
          surface_primary = NA_real_,
          surface_rmse = NA_real_,
          surface_mae = NA_real_,
          surface_iae = NA_real_,
          surface_ise = NA_real_,
          ok = FALSE,
          error = conditionMessage(e),
          stringsAsFactors = FALSE
        ),
        eval = NULL
      )
    }
  )
  res
}

parse_bandwidth_string <- function(x) {
  if (is.na(x) || !nzchar(x)) return(numeric(0))
  as.numeric(strsplit(x, ";", fixed = TRUE)[[1L]])
}

bandwidth_l2_diff <- function(a, b) {
  va <- parse_bandwidth_string(a)
  vb <- parse_bandwidth_string(b)
  if (length(va) == 0L || length(vb) == 0L || length(va) != length(vb)) return(NA_real_)
  sqrt(sum((va - vb)^2))
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  mean(x)
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

summarize_results <- function(raw_runs, raw_eval, nmulti_ref) {
  ok_runs <- raw_runs[raw_runs$ok, , drop = FALSE]
  if (nrow(ok_runs) == 0L) return(list(summary = data.frame(), drift = data.frame()))

  key_seed <- interaction(ok_runs$scenario_id, ok_runs$n, ok_runs$seed, drop = TRUE)
  best_primary <- tapply(ok_runs$surface_primary, key_seed, min, na.rm = TRUE)
  ok_runs$best_primary <- as.numeric(best_primary[key_seed])
  ok_runs$within_1pct_best <- ok_runs$surface_primary <= 1.01 * ok_runs$best_primary

  ref_runs <- ok_runs[ok_runs$nmulti == nmulti_ref, , drop = FALSE]
  names(ref_runs)[names(ref_runs) %in% c("elapsed_bw", "elapsed_total", "fval", "bandwidth", "surface_primary")] <-
    c("elapsed_bw_ref", "elapsed_total_ref", "fval_ref", "bandwidth_ref", "surface_primary_ref")

  merged <- merge(
    ok_runs,
    ref_runs[, c("scenario_id", "n", "seed", "elapsed_bw_ref", "elapsed_total_ref", "fval_ref", "bandwidth_ref", "surface_primary_ref")],
    by = c("scenario_id", "n", "seed"),
    all.x = TRUE,
    sort = FALSE
  )
  merged$bw_speedup_vs_ref <- merged$elapsed_bw_ref / merged$elapsed_bw
  merged$total_speedup_vs_ref <- merged$elapsed_total_ref / merged$elapsed_total
  merged$fval_abs_gap_vs_ref <- abs(merged$fval - merged$fval_ref)
  merged$bandwidth_l2_vs_ref <- mapply(bandwidth_l2_diff, merged$bandwidth, merged$bandwidth_ref)
  merged$surface_primary_delta_vs_ref <- merged$surface_primary - merged$surface_primary_ref

  ref_eval <- raw_eval[raw_eval$nmulti == nmulti_ref, , drop = FALSE]
  names(ref_eval)[names(ref_eval) == "pred"] <- "pred_ref"
  eval_merged <- merge(
    raw_eval,
    ref_eval[, c("scenario_id", "n", "seed", "eval_id", "pred_ref")],
    by = c("scenario_id", "n", "seed", "eval_id"),
    all.x = TRUE,
    sort = FALSE
  )

  drift_groups <- split(eval_merged, interaction(eval_merged$scenario_id, eval_merged$n, eval_merged$seed, eval_merged$nmulti, drop = TRUE))
  drift_rows <- lapply(drift_groups, function(df) {
    data.frame(
      scenario_id = df$scenario_id[1L],
      family = df$family[1L],
      regtype = df$regtype[1L],
      degree = df$degree[1L],
      n = df$n[1L],
      seed = df$seed[1L],
      nmulti = df$nmulti[1L],
      surface_drift_rmse_vs_ref = sqrt(mean((df$pred - df$pred_ref)^2)),
      surface_drift_max_abs_vs_ref = max(abs(df$pred - df$pred_ref)),
      stringsAsFactors = FALSE
    )
  })
  drift <- do.call(rbind, drift_rows)

  merged <- merge(
    merged,
    drift[, c("scenario_id", "n", "seed", "nmulti", "surface_drift_rmse_vs_ref", "surface_drift_max_abs_vs_ref")],
    by = c("scenario_id", "n", "seed", "nmulti"),
    all.x = TRUE,
    sort = FALSE
  )

  group_key <- interaction(merged$scenario_id, merged$n, merged$nmulti, drop = TRUE)
  summary_rows <- lapply(split(seq_len(nrow(merged)), group_key), function(idx) {
    df <- merged[idx, , drop = FALSE]
    data.frame(
      scenario_id = df$scenario_id[1L],
      family = df$family[1L],
      regtype = df$regtype[1L],
      degree = df$degree[1L],
      n = df$n[1L],
      nmulti = df$nmulti[1L],
      runs = nrow(df),
      median_elapsed_bw = safe_median(df$elapsed_bw),
      median_elapsed_total = safe_median(df$elapsed_total),
      median_bw_speedup_vs_ref = safe_median(df$bw_speedup_vs_ref),
      median_total_speedup_vs_ref = safe_median(df$total_speedup_vs_ref),
      median_surface_primary = safe_median(df$surface_primary),
      median_surface_primary_delta_vs_ref = safe_median(df$surface_primary_delta_vs_ref),
      median_surface_drift_rmse_vs_ref = safe_median(df$surface_drift_rmse_vs_ref),
      p90_surface_drift_rmse_vs_ref = safe_quantile(df$surface_drift_rmse_vs_ref, 0.90),
      median_surface_drift_max_abs_vs_ref = safe_median(df$surface_drift_max_abs_vs_ref),
      median_fval_abs_gap_vs_ref = safe_median(df$fval_abs_gap_vs_ref),
      median_bandwidth_l2_vs_ref = safe_median(df$bandwidth_l2_vs_ref),
      share_within_1pct_best = safe_mean(df$within_1pct_best),
      stringsAsFactors = FALSE
    )
  })

  list(summary = do.call(rbind, summary_rows), drift = drift)
}

write_report <- function(summary_df, out_file, nmulti_ref) {
  lines <- c(
    "# nmulti Core Surface Study",
    "",
    sprintf("Reference `nmulti`: `%d`.", nmulti_ref),
    "",
    "Decision lens:",
    "- lower `nmulti` is attractive only if runtime savings are meaningful and surface drift versus the reference is tiny;",
    "- `share_within_1pct_best` measures how often a candidate is within 1% of the best seed-matched truth-based surface error;",
    "- raw bandwidth drift is diagnostic, not decisive, because very different bandwidths can still induce nearly identical fitted surfaces.",
    ""
  )

  families <- split(summary_df, summary_df$family)
  for (fam in names(families)) {
    lines <- c(lines, sprintf("## %s", fam), "")
    by_scenario <- split(families[[fam]], families[[fam]]$scenario_id)
    for (sid in names(by_scenario)) {
      df <- by_scenario[[sid]]
      df <- df[order(df$n, df$nmulti), , drop = FALSE]
      lines <- c(lines, sprintf("### %s", sid), "")
      lines <- c(lines, "| n | nmulti | med bw sec | med bw speedup | med surf drift RMSE | share within 1% best |")
      lines <- c(lines, "|---:|---:|---:|---:|---:|---:|")
      for (i in seq_len(nrow(df))) {
        lines <- c(lines, sprintf(
          "| %d | %d | %.3f | %.2f | %.6f | %.2f |",
          df$n[i],
          df$nmulti[i],
          df$median_elapsed_bw[i],
          df$median_bw_speedup_vs_ref[i],
          df$median_surface_drift_rmse_vs_ref[i],
          df$share_within_1pct_best[i]
        ))
      }
      lines <- c(lines, "")
    }
  }

  writeLines(lines, con = out_file)
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  cfg <- parse_args(args)
  dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
  load_backend()

  scenarios <- scenario_catalog()
  reg_eval <- build_regression_eval()
  dens_eval <- build_density_eval()

  raw_runs <- list()
  raw_eval <- list()
  idx <- 0L
  total_jobs <- length(scenarios) * length(cfg$n_grid) * length(cfg$seeds) * length(cfg$nmulti_grid)

  if (cfg$show_progress) {
    cat("Running nmulti core surface study\n")
    cat("out_dir=", cfg$out_dir, "\n", sep = "")
    cat("n_grid=", paste(cfg$n_grid, collapse = ","), "\n", sep = "")
    cat("seeds=", paste(cfg$seeds, collapse = ","), "\n", sep = "")
    cat("nmulti_grid=", paste(cfg$nmulti_grid, collapse = ","), "\n", sep = "")
  }

  for (spec in scenarios) {
    eval_df <- if (identical(spec$family, "npreg")) reg_eval else dens_eval
    for (n in cfg$n_grid) {
      for (seed in cfg$seeds) {
        for (nmulti in cfg$nmulti_grid) {
          idx <- idx + 1L
          if (cfg$show_progress) {
            cat(sprintf("[%d/%d] %s n=%d seed=%d nmulti=%d\n", idx, total_jobs, spec$id, n, seed, nmulti))
          }
          res <- run_one(spec = spec, n = n, seed = seed, nmulti = nmulti, eval_df = eval_df)
          raw_runs[[length(raw_runs) + 1L]] <- res$run
          if (!is.null(res$eval)) raw_eval[[length(raw_eval) + 1L]] <- res$eval
        }
      }
    }
  }

  raw_runs_df <- do.call(rbind, raw_runs)
  raw_eval_df <- do.call(rbind, raw_eval)
  nmulti_ref <- max(cfg$nmulti_grid)
  summaries <- summarize_results(raw_runs_df, raw_eval_df, nmulti_ref = nmulti_ref)

  runs_file <- file.path(cfg$out_dir, "nmulti_core_surface_runs.csv")
  eval_file <- file.path(cfg$out_dir, "nmulti_core_surface_eval.csv")
  summary_file <- file.path(cfg$out_dir, "nmulti_core_surface_summary.csv")
  drift_file <- file.path(cfg$out_dir, "nmulti_core_surface_drift.csv")
  report_file <- file.path(cfg$out_dir, "nmulti_core_surface_report.md")

  write.csv(raw_runs_df, runs_file, row.names = FALSE)
  write.csv(raw_eval_df, eval_file, row.names = FALSE)
  write.csv(summaries$summary, summary_file, row.names = FALSE)
  write.csv(summaries$drift, drift_file, row.names = FALSE)
  write_report(summaries$summary, out_file = report_file, nmulti_ref = nmulti_ref)

  if (cfg$show_progress) {
    cat("Runs:", runs_file, "\n")
    cat("Eval:", eval_file, "\n")
    cat("Summary:", summary_file, "\n")
    cat("Drift:", drift_file, "\n")
    cat("Report:", report_file, "\n")
    print(summaries$summary)
  }
}

if (sys.nframe() == 0L) main()
