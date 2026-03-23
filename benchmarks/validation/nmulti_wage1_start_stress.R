#!/usr/bin/env Rscript

parse_args <- function(args) {
  out <- list(
    degrees = c(1L, 2L),
    nmulti_grid = c(1L, 2L, 5L),
    cat_starts = c(1e-6, 0.05, 0.25, 0.50, 0.90),
    cont_multipliers = c(0.25, 0.50, 1.0, 2.0),
    out_dir = sprintf("/tmp/nmulti_wage1_start_stress_%s", format(Sys.time(), "%Y%m%d_%H%M%S")),
    show_progress = TRUE
  )

  parse_int_csv <- function(x, name) {
    vals <- as.integer(strsplit(x, ",", fixed = TRUE)[[1L]])
    if (length(vals) == 0L || anyNA(vals)) stop(name, " must be a comma-separated integer list")
    vals
  }

  parse_num_csv <- function(x, name) {
    vals <- as.numeric(strsplit(x, ",", fixed = TRUE)[[1L]])
    if (length(vals) == 0L || anyNA(vals)) stop(name, " must be a comma-separated numeric list")
    vals
  }

  if (length(args) == 0L) return(out)

  for (a in args) {
    if (!grepl("^--", a)) stop("Bad arg: ", a)
    raw <- sub("^--", "", a)
    key <- sub("=.*$", "", raw)
    val <- if (grepl("=", raw, fixed = TRUE)) sub("^[^=]*=", "", raw) else "TRUE"

    if (key == "degrees") out$degrees <- parse_int_csv(val, "degrees")
    else if (key == "nmulti_grid") out$nmulti_grid <- parse_int_csv(val, "nmulti_grid")
    else if (key == "cat_starts") out$cat_starts <- parse_num_csv(val, "cat_starts")
    else if (key == "cont_multipliers") out$cont_multipliers <- parse_num_csv(val, "cont_multipliers")
    else if (key == "out_dir") out$out_dir <- val
    else if (key == "show_progress") out$show_progress <- as.logical(val)
    else stop("Unknown arg: ", key)
  }

  out$degrees <- sort(unique(out$degrees))
  out$nmulti_grid <- sort(unique(out$nmulti_grid))
  out$cat_starts <- sort(unique(out$cat_starts))
  out$cont_multipliers <- sort(unique(out$cont_multipliers))

  if (any(!out$degrees %in% c(1L, 2L))) stop("degrees must be drawn from 1,2")
  if (length(out$nmulti_grid) < 2L || !5L %in% out$nmulti_grid) stop("nmulti_grid must include 5 as the reference")
  if (any(out$nmulti_grid < 1L)) stop("nmulti_grid values must be >= 1")
  if (any(out$cat_starts <= 0 | out$cat_starts >= 1)) stop("cat_starts must lie strictly between 0 and 1")
  if (any(out$cont_multipliers <= 0)) stop("cont_multipliers must be positive")

  out
}

load_backend <- function() {
  suppressPackageStartupMessages(library(np))
  data("wage1", package = "np")
  options(np.messages = FALSE)
}

fit_formula <- function() {
  lwage ~ married + female + nonwhite + educ + exper + tenure
}

fit_reference <- function(degree) {
  set.seed(20260322 + degree)
  npregbw(
    formula = fit_formula(),
    regtype = "lp",
    degree = rep.int(degree, 3L),
    data = wage1,
    nmulti = 5L
  )
}

make_start_grid <- function(ref_bw, cat_starts, cont_multipliers) {
  ref_vals <- as.numeric(ref_bw$bw)
  idx_cat <- seq_len(3L)
  idx_con <- 4:6
  starts <- vector("list", length(cat_starts) * length(cont_multipliers))
  k <- 0L

  for (cat_start in cat_starts) {
    for (cont_mult in cont_multipliers) {
      k <- k + 1L
      cur <- ref_vals
      cur[idx_cat] <- rep.int(cat_start, length(idx_cat))
      cur[idx_con] <- pmax(ref_vals[idx_con] * cont_mult, 1e-8)
      starts[[k]] <- list(
        start_id = sprintf("cat_%s__cont_%s",
                           format(signif(cat_start, 3L), scientific = FALSE, trim = TRUE),
                           format(signif(cont_mult, 3L), scientific = FALSE, trim = TRUE)),
        cat_start = cat_start,
        cont_multiplier = cont_mult,
        bws = cur
      )
    }
  }

  starts
}

encode_numeric <- function(x) paste(signif(as.numeric(x), 10L), collapse = ";")

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

fitted_drift <- function(a, b) {
  list(
    rmse = sqrt(mean((a - b)^2)),
    max_abs = max(abs(a - b))
  )
}

run_one <- function(degree, start_spec, nmulti, seed_base) {
  set.seed(seed_base + nmulti)
  t0 <- proc.time()[["elapsed"]]
  bw <- npregbw(
    formula = fit_formula(),
    regtype = "lp",
    degree = rep.int(degree, 3L),
    data = wage1,
    bws = start_spec$bws,
    nmulti = nmulti
  )
  elapsed_bw <- proc.time()[["elapsed"]] - t0

  fit <- npreg(bws = bw, data = wage1)

  list(
    degree = degree,
    start_id = start_spec$start_id,
    cat_start = start_spec$cat_start,
    cont_multiplier = start_spec$cont_multiplier,
    nmulti = nmulti,
    elapsed_bw = elapsed_bw,
    fval = first_num_or_na(bw$fval),
    ifval = first_num_or_na(bw$ifval),
    num_fval = first_num_or_na(bw$num.fval),
    eval_count = count_history(bw$eval.history),
    invalid_count = count_history(bw$invalid.history),
    bandwidth = encode_numeric(bw$bw),
    fitted = fitted(fit)
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

  raw_rows <- list()
  drift_rows <- list()
  winner_rows <- list()
  raw_idx <- 0L
  drift_idx <- 0L
  winner_idx <- 0L

  total_jobs <- length(cfg$degrees) * length(cfg$cat_starts) * length(cfg$cont_multipliers) * length(cfg$nmulti_grid)
  job_idx <- 0L

  if (cfg$show_progress) {
    cat("Running wage1 start-value stress test\n")
    cat("out_dir=", cfg$out_dir, "\n", sep = "")
    cat("degrees=", paste(cfg$degrees, collapse = ","), "\n", sep = "")
    cat("nmulti_grid=", paste(cfg$nmulti_grid, collapse = ","), "\n", sep = "")
    cat("cat_starts=", paste(cfg$cat_starts, collapse = ","), "\n", sep = "")
    cat("cont_multipliers=", paste(cfg$cont_multipliers, collapse = ","), "\n", sep = "")
  }

  for (degree in cfg$degrees) {
    ref_bw <- fit_reference(degree)
    starts <- make_start_grid(ref_bw, cfg$cat_starts, cfg$cont_multipliers)

    for (i in seq_along(starts)) {
      start_spec <- starts[[i]]
      seed_base <- 7000L + 100L * degree + i
      per_start <- vector("list", length(cfg$nmulti_grid))

      for (j in seq_along(cfg$nmulti_grid)) {
        nmulti <- cfg$nmulti_grid[j]
        job_idx <- job_idx + 1L
        if (cfg$show_progress) {
          cat(sprintf("[%d/%d] degree=%d start=%s nmulti=%d\n",
                      job_idx, total_jobs, degree, start_spec$start_id, nmulti))
        }
        per_start[[j]] <- run_one(degree = degree, start_spec = start_spec, nmulti = nmulti, seed_base = seed_base)
      }

      ref <- per_start[[which(cfg$nmulti_grid == 5L)]]
      fvals <- vapply(per_start, `[[`, numeric(1L), "fval")
      winner_nm <- per_start[[which.min(fvals)]]$nmulti

      for (res in per_start) {
        raw_idx <- raw_idx + 1L
        d <- fitted_drift(res$fitted, ref$fitted)
        raw_rows[[raw_idx]] <- data.frame(
          degree = res$degree,
          start_id = res$start_id,
          cat_start = res$cat_start,
          cont_multiplier = res$cont_multiplier,
          nmulti = res$nmulti,
          elapsed_bw = res$elapsed_bw,
          fval = res$fval,
          ifval = res$ifval,
          num_fval = res$num_fval,
          eval_count = res$eval_count,
          invalid_count = res$invalid_count,
          bandwidth = res$bandwidth,
          bandwidth_l2_vs_ref5 = bandwidth_l2_diff(res$bandwidth, ref$bandwidth),
          fval_gap_vs_ref5 = abs(res$fval - ref$fval),
          fitted_rmse_vs_ref5 = d$rmse,
          fitted_max_abs_vs_ref5 = d$max_abs,
          winner_nm = winner_nm,
          within_1pct_best_fval = res$fval <= 1.01 * min(fvals),
          stringsAsFactors = FALSE
        )
      }

      winner_idx <- winner_idx + 1L
      winner_rows[[winner_idx]] <- data.frame(
        degree = degree,
        start_id = start_spec$start_id,
        cat_start = start_spec$cat_start,
        cont_multiplier = start_spec$cont_multiplier,
        best_nm = winner_nm,
        best_fval = min(fvals),
        stringsAsFactors = FALSE
      )
    }
  }

  raw_df <- do.call(rbind, raw_rows)
  winner_df <- do.call(rbind, winner_rows)

  split_key <- interaction(raw_df$degree, raw_df$nmulti, drop = TRUE)
  summary_rows <- lapply(split(seq_len(nrow(raw_df)), split_key), function(idx) {
    df <- raw_df[idx, , drop = FALSE]
    data.frame(
      degree = df$degree[1L],
      nmulti = df$nmulti[1L],
      runs = nrow(df),
      median_elapsed_bw = safe_median(df$elapsed_bw),
      median_bw_speedup_vs_ref5 = safe_median(raw_df$elapsed_bw[raw_df$degree == df$degree[1L] & raw_df$nmulti == 5L]) / safe_median(df$elapsed_bw),
      median_fval_gap_vs_ref5 = safe_median(df$fval_gap_vs_ref5),
      p90_fval_gap_vs_ref5 = safe_quantile(df$fval_gap_vs_ref5, 0.90),
      median_fitted_rmse_vs_ref5 = safe_median(df$fitted_rmse_vs_ref5),
      p90_fitted_rmse_vs_ref5 = safe_quantile(df$fitted_rmse_vs_ref5, 0.90),
      median_fitted_max_abs_vs_ref5 = safe_median(df$fitted_max_abs_vs_ref5),
      p90_fitted_max_abs_vs_ref5 = safe_quantile(df$fitted_max_abs_vs_ref5, 0.90),
      median_bandwidth_l2_vs_ref5 = safe_median(df$bandwidth_l2_vs_ref5),
      share_best_fval = mean(df$winner_nm == df$nmulti[1L]),
      share_within_1pct_best_fval = mean(df$within_1pct_best_fval),
      stringsAsFactors = FALSE
    )
  })
  summary_df <- do.call(rbind, summary_rows)

  runs_file <- file.path(cfg$out_dir, "nmulti_wage1_start_stress_runs.csv")
  winners_file <- file.path(cfg$out_dir, "nmulti_wage1_start_stress_winners.csv")
  summary_file <- file.path(cfg$out_dir, "nmulti_wage1_start_stress_summary.csv")
  report_file <- file.path(cfg$out_dir, "nmulti_wage1_start_stress_report.md")

  write.csv(raw_df, runs_file, row.names = FALSE)
  write.csv(winner_df, winners_file, row.names = FALSE)
  write.csv(summary_df, summary_file, row.names = FALSE)

  report_lines <- c(
    "# wage1 LP Start-Value Stress Test",
    "",
    "This study fixes the `wage1` LP regression route and varies the starting bandwidth vector.",
    "Each start cell is run with `nmulti = 1, 2, 5`.",
    "`nmulti = 5` is the practical reference, and lower `nmulti` values are judged by:",
    "- objective gap relative to the same-start `nmulti=5` run,",
    "- fitted-value drift relative to the same-start `nmulti=5` run,",
    "- how often they achieve the best `fval` across `1, 2, 5`.",
    ""
  )

  for (deg in cfg$degrees) {
    sub <- summary_df[summary_df$degree == deg, , drop = FALSE]
    sub <- sub[order(sub$nmulti), , drop = FALSE]
    report_lines <- c(report_lines, sprintf("## degree=%d", deg), "")
    report_lines <- c(report_lines, "| nmulti | med bw sec | med speedup vs 5 | med fval gap | p90 fitted RMSE vs 5 | share best fval | share within 1% best fval |")
    report_lines <- c(report_lines, "|---:|---:|---:|---:|---:|---:|---:|")
    for (i in seq_len(nrow(sub))) {
      report_lines <- c(report_lines, sprintf(
        "| %d | %.3f | %.2f | %.6g | %.6g | %.2f | %.2f |",
        sub$nmulti[i],
        sub$median_elapsed_bw[i],
        sub$median_bw_speedup_vs_ref5[i],
        sub$median_fval_gap_vs_ref5[i],
        sub$p90_fitted_rmse_vs_ref5[i],
        sub$share_best_fval[i],
        sub$share_within_1pct_best_fval[i]
      ))
    }
    report_lines <- c(report_lines, "")
  }
  writeLines(report_lines, con = report_file)

  if (cfg$show_progress) {
    cat("Runs:", runs_file, "\n")
    cat("Winners:", winners_file, "\n")
    cat("Summary:", summary_file, "\n")
    cat("Report:", report_file, "\n")
    print(summary_df)
  }
}

if (sys.nframe() == 0L) main()
