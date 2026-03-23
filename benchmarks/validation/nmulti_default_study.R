#!/usr/bin/env Rscript

parse_args <- function(args) {
  out <- list(
    preset = "smoke",
    scenarios = NULL,
    n_grid = NULL,
    nmulti_grid = c(1L, 2L, 5L),
    seeds = NULL,
    out_dir = sprintf("/tmp/nmulti_default_study_%s", format(Sys.time(), "%Y%m%d_%H%M%S")),
    show_progress = TRUE
  )

  parse_int_csv <- function(x, name) {
    vals <- as.integer(strsplit(x, ",", fixed = TRUE)[[1L]])
    if (length(vals) == 0L || anyNA(vals)) {
      stop(name, " must be a comma-separated list of integers")
    }
    vals
  }

  parse_chr_csv <- function(x) {
    vals <- strsplit(x, ",", fixed = TRUE)[[1L]]
    vals[nzchar(vals)]
  }

  if (length(args) == 0L) return(out)

  for (a in args) {
    if (!grepl("^--", a)) stop("Bad arg: ", a)
    raw <- sub("^--", "", a)
    key <- sub("=.*$", "", raw)
    val <- if (grepl("=", raw, fixed = TRUE)) sub("^[^=]*=", "", raw) else "TRUE"

    if (key == "preset") out$preset <- val
    else if (key == "scenarios") out$scenarios <- parse_chr_csv(val)
    else if (key == "n_grid") out$n_grid <- parse_int_csv(val, "n_grid")
    else if (key == "nmulti_grid") out$nmulti_grid <- parse_int_csv(val, "nmulti_grid")
    else if (key == "seeds") out$seeds <- parse_int_csv(val, "seeds")
    else if (key == "out_dir") out$out_dir <- val
    else if (key == "show_progress") out$show_progress <- as.logical(val)
    else stop("Unknown arg: ", key)
  }

  if (!out$preset %in% c("smoke", "full")) {
    stop("preset must be smoke or full")
  }
  if (length(out$nmulti_grid) < 2L) {
    stop("nmulti_grid must contain at least two values")
  }
  if (any(out$nmulti_grid < 1L)) {
    stop("nmulti_grid values must be >= 1")
  }

  out$nmulti_grid <- sort(unique(out$nmulti_grid))
  out
}

load_backend <- function() {
  suppressPackageStartupMessages(library(np))
  options(np.messages = FALSE)
}

scenario_catalog <- function() {
  list(
    npreg_ll_fixed_mix = list(
      id = "npreg_ll_fixed_mix",
      family = "npreg",
      predictors = c("x1", "x2", "z1", "z2"),
      regtype = "ll",
      bwtype = "fixed",
      bwmethod = "cv.ls"
    ),
    npreg_ll_generalized_nn_mix = list(
      id = "npreg_ll_generalized_nn_mix",
      family = "npreg",
      predictors = c("x1", "x2", "z1", "z2"),
      regtype = "ll",
      bwtype = "generalized_nn",
      bwmethod = "cv.ls"
    ),
    npreg_lp_fixed_mix = list(
      id = "npreg_lp_fixed_mix",
      family = "npreg",
      predictors = c("x1", "x2", "z1", "z2"),
      regtype = "lp",
      bwtype = "fixed",
      bwmethod = "cv.ls",
      basis = "glp",
      degree = c(2L, 2L),
      bernstein.basis = FALSE
    ),
    npcdens_lc_fixed_mix = list(
      id = "npcdens_lc_fixed_mix",
      family = "npcdens",
      predictors = c("x1", "x2", "z1"),
      regtype = "lc",
      bwtype = "fixed",
      bwmethod = "cv.ml"
    ),
    npcdens_lc_generalized_nn_mix = list(
      id = "npcdens_lc_generalized_nn_mix",
      family = "npcdens",
      predictors = c("x1", "x2", "z1"),
      regtype = "lc",
      bwtype = "generalized_nn",
      bwmethod = "cv.ml"
    ),
    npcdens_lp_fixed_mix = list(
      id = "npcdens_lp_fixed_mix",
      family = "npcdens",
      predictors = c("x1", "x2", "z1"),
      regtype = "lp",
      bwtype = "fixed",
      bwmethod = "cv.ml",
      basis = "glp",
      degree = c(2L, 2L),
      bernstein.basis = FALSE
    ),
    npcdist_lc_fixed_mix = list(
      id = "npcdist_lc_fixed_mix",
      family = "npcdist",
      predictors = c("x1", "x2", "z1"),
      regtype = "lc",
      bwtype = "fixed",
      bwmethod = "cv.ls"
    ),
    npcdist_lc_generalized_nn_mix = list(
      id = "npcdist_lc_generalized_nn_mix",
      family = "npcdist",
      predictors = c("x1", "x2", "z1"),
      regtype = "lc",
      bwtype = "generalized_nn",
      bwmethod = "cv.ls"
    ),
    npcdist_lp_fixed_mix = list(
      id = "npcdist_lp_fixed_mix",
      family = "npcdist",
      predictors = c("x1", "x2", "z1"),
      regtype = "lp",
      bwtype = "fixed",
      bwmethod = "cv.ls",
      basis = "glp",
      degree = c(2L, 2L),
      bernstein.basis = FALSE
    )
  )
}

resolve_config <- function(cfg) {
  catalog <- scenario_catalog()
  smoke_ids <- c("npreg_ll_fixed_mix", "npcdens_lc_fixed_mix", "npcdist_lc_fixed_mix")
  full_ids <- names(catalog)

  scenario_ids <- if (is.null(cfg$scenarios)) {
    if (cfg$preset == "smoke") smoke_ids else full_ids
  } else {
    cfg$scenarios
  }

  missing_ids <- setdiff(scenario_ids, names(catalog))
  if (length(missing_ids) > 0L) {
    stop("Unknown scenario ids: ", paste(missing_ids, collapse = ", "))
  }

  n_grid <- if (!is.null(cfg$n_grid)) {
    cfg$n_grid
  } else if (cfg$preset == "smoke") {
    c(100L, 250L)
  } else {
    c(100L, 250L, 500L)
  }

  seeds <- if (!is.null(cfg$seeds)) {
    cfg$seeds
  } else if (cfg$preset == "smoke") {
    101:103
  } else {
    101:120
  }

  list(
    scenarios = unname(catalog[scenario_ids]),
    n_grid = sort(unique(as.integer(n_grid))),
    nmulti_grid = sort(unique(as.integer(cfg$nmulti_grid))),
    seeds = sort(unique(as.integer(seeds))),
    out_dir = cfg$out_dir,
    show_progress = isTRUE(cfg$show_progress)
  )
}

generate_data <- function(seed, n) {
  set.seed(seed)
  x1 <- runif(n)
  x2 <- runif(n)
  z1 <- factor(sample(c("a", "b", "c"), n, replace = TRUE, prob = c(0.35, 0.40, 0.25)),
               levels = c("a", "b", "c"))
  z2 <- ordered(sample(c("1", "2", "3", "4"), n, replace = TRUE, prob = c(0.20, 0.25, 0.30, 0.25)),
                levels = c("1", "2", "3", "4"))

  mu <- truth_mu(data.frame(x1 = x1, x2 = x2, z1 = z1, z2 = z2))
  sigma <- truth_sigma(data.frame(x1 = x1, x2 = x2, z1 = z1, z2 = z2))
  y <- mu + sigma * rnorm(n)

  data.frame(y = y, x1 = x1, x2 = x2, z1 = z1, z2 = z2)
}

truth_mu <- function(df) {
  n <- nrow(df)
  x1 <- if ("x1" %in% names(df)) df$x1 else rep(0.5, n)
  x2 <- if ("x2" %in% names(df)) df$x2 else rep(0.5, n)
  z1 <- if ("z1" %in% names(df)) as.character(df$z1) else rep("a", n)
  z2 <- if ("z2" %in% names(df)) as.character(df$z2) else rep("1", n)

  z1_eff <- c(a = 0.00, b = 0.45, c = -0.30)[z1]
  z2_eff <- c("1" = 0.00, "2" = 0.20, "3" = 0.45, "4" = 0.75)[z2]
  0.75 * sin(2 * pi * x1) + 0.50 * cos(pi * x2) + z1_eff + z2_eff
}

truth_sigma <- function(df) {
  n <- nrow(df)
  x2 <- if ("x2" %in% names(df)) df$x2 else rep(0.5, n)
  0.25 + 0.15 * x2
}

make_formula <- function(spec) {
  as.formula(paste("y ~", paste(spec$predictors, collapse = " + ")))
}

subset_predictors <- function(df, predictors) {
  keep <- intersect(c("y", predictors), names(df))
  df[, keep, drop = FALSE]
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

build_x_anchors <- function(spec) {
  anchors <- data.frame(
    x1 = c(0.15, 0.20, 0.50, 0.85, 0.80, 0.35),
    x2 = c(0.20, 0.85, 0.50, 0.15, 0.80, 0.35),
    z1 = factor(c("a", "b", "c", "a", "c", "b"), levels = c("a", "b", "c")),
    z2 = ordered(c("1", "2", "3", "4", "2", "3"), levels = c("1", "2", "3", "4")),
    stringsAsFactors = FALSE
  )
  anchors[, spec$predictors, drop = FALSE]
}

build_eval_data <- function(spec) {
  anchors <- build_x_anchors(spec)

  if (identical(spec$family, "npreg")) {
    out <- anchors
    out$eval_id <- seq_len(nrow(out))
    return(out)
  }

  y_grid <- seq(-3.5, 4.5, length.out = 81L)
  rows <- vector("list", nrow(anchors))
  for (i in seq_len(nrow(anchors))) {
    row_i <- anchors[rep(i, length(y_grid)), , drop = FALSE]
    row_i$y <- y_grid
    row_i$slice_id <- i
    row_i$eval_id <- seq_len(length(y_grid)) + (i - 1L) * length(y_grid)
    rows[[i]] <- row_i
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[, c("y", spec$predictors, "slice_id", "eval_id"), drop = FALSE]
}

trapz <- function(x, y) {
  n <- length(x)
  if (n < 2L) return(NA_real_)
  sum((x[-1L] - x[-n]) * (y[-1L] + y[-n]) / 2)
}

evaluate_truth <- function(spec, eval_data) {
  if (identical(spec$family, "npreg")) {
    return(truth_mu(eval_data))
  }

  x_only <- eval_data[, spec$predictors, drop = FALSE]
  mu <- truth_mu(x_only)
  sigma <- truth_sigma(x_only)
  y <- eval_data$y

  if (identical(spec$family, "npcdens")) {
    return(dnorm(y, mean = mu, sd = sigma))
  }

  pnorm(y, mean = mu, sd = sigma)
}

compute_eval_metrics <- function(spec, eval_data, pred, truth) {
  err <- pred - truth
  out <- list(
    pred_rmse = sqrt(mean(err^2)),
    pred_mae = mean(abs(err)),
    pred_max_abs = max(abs(err))
  )

  if (!identical(spec$family, "npreg")) {
    slices <- split(seq_len(nrow(eval_data)), eval_data$slice_id)
    iae <- numeric(length(slices))
    ise <- numeric(length(slices))
    for (i in seq_along(slices)) {
      idx <- slices[[i]]
      y <- eval_data$y[idx]
      iae[i] <- trapz(y, abs(err[idx]))
      ise[i] <- trapz(y, err[idx]^2)
    }
    out$grid_iae <- mean(iae)
    out$grid_ise <- mean(ise)
  } else {
    out$grid_iae <- NA_real_
    out$grid_ise <- NA_real_
  }

  out
}

fit_case <- function(spec, dat, nmulti, eval_template) {
  bw_fun <- switch(
    spec$family,
    npreg = npregbw,
    npcdens = npcdensbw,
    npcdist = npcdistbw
  )
  fit_fun <- switch(
    spec$family,
    npreg = npreg,
    npcdens = npcdens,
    npcdist = npcdist
  )

  fit_args <- list(
    formula = make_formula(spec),
    data = dat,
    nmulti = nmulti,
    bwtype = spec$bwtype,
    bwmethod = spec$bwmethod,
    regtype = spec$regtype
  )

  if (!is.null(spec$basis)) fit_args$basis <- spec$basis
  if (!is.null(spec$degree)) fit_args$degree <- spec$degree
  if (!is.null(spec$bernstein.basis)) fit_args$bernstein.basis <- spec$bernstein.basis

  t0 <- proc.time()[["elapsed"]]
  bw <- do.call(bw_fun, fit_args)
  elapsed_bw <- proc.time()[["elapsed"]] - t0

  t1 <- proc.time()[["elapsed"]]
  fit <- fit_fun(bws = bw, data = dat)
  elapsed_fit <- proc.time()[["elapsed"]] - t1

  eval_newdata <- eval_template
  pred_data <- eval_newdata[, intersect(c("y", spec$predictors), names(eval_newdata)), drop = FALSE]
  pred <- as.numeric(predict(fit, newdata = pred_data))
  truth <- evaluate_truth(spec, eval_newdata)
  metrics <- compute_eval_metrics(spec, eval_newdata, pred, truth)

  eval_df <- data.frame(
    scenario_id = spec$id,
    family = spec$family,
    n = nrow(dat),
    seed = NA_integer_,
    nmulti = nmulti,
    eval_id = eval_newdata$eval_id,
    slice_id = if ("slice_id" %in% names(eval_newdata)) eval_newdata$slice_id else NA_integer_,
    pred = pred,
    truth = truth,
    stringsAsFactors = FALSE
  )

  list(
    bw = bw,
    fit = fit,
    elapsed_bw = elapsed_bw,
    elapsed_fit = elapsed_fit,
    metrics = metrics,
    eval_df = eval_df
  )
}

run_one <- function(spec, n, seed, nmulti, eval_template) {
  dat <- generate_data(seed = seed, n = n)

  out <- tryCatch(
    {
      res <- fit_case(spec = spec, dat = dat, nmulti = nmulti, eval_template = eval_template)
      bw_vals <- extract_bw_values(res$bw)
      run_row <- data.frame(
        scenario_id = spec$id,
        family = spec$family,
        n = n,
        seed = seed,
        nmulti = nmulti,
        regtype = spec$regtype,
        bwtype = spec$bwtype,
        bwmethod = spec$bwmethod,
        elapsed_bw = res$elapsed_bw,
        elapsed_fit = res$elapsed_fit,
        elapsed_total = res$elapsed_bw + res$elapsed_fit,
        fval = first_num_or_na(res$bw$fval),
        ifval = first_num_or_na(res$bw$ifval),
        eval_count = count_history(res$bw$eval.history),
        invalid_count = count_history(res$bw$invalid.history),
        num_fval = first_num_or_na(if (!is.null(res$bw$num.fval)) res$bw$num.fval else res$bw$num.feval),
        bandwidth = encode_numeric(bw_vals),
        pred_rmse = res$metrics$pred_rmse,
        pred_mae = res$metrics$pred_mae,
        pred_max_abs = res$metrics$pred_max_abs,
        grid_iae = res$metrics$grid_iae,
        grid_ise = res$metrics$grid_ise,
        ok = TRUE,
        error = "",
        stringsAsFactors = FALSE
      )

      res$eval_df$seed <- seed
      list(run = run_row, eval = res$eval_df)
    },
    error = function(e) {
      list(
        run = data.frame(
          scenario_id = spec$id,
          family = spec$family,
          n = n,
          seed = seed,
          nmulti = nmulti,
          regtype = spec$regtype,
          bwtype = spec$bwtype,
          bwmethod = spec$bwmethod,
          elapsed_bw = NA_real_,
          elapsed_fit = NA_real_,
          elapsed_total = NA_real_,
          fval = NA_real_,
          ifval = NA_real_,
          eval_count = NA_integer_,
          invalid_count = NA_integer_,
          num_fval = NA_real_,
          bandwidth = "",
          pred_rmse = NA_real_,
          pred_mae = NA_real_,
          pred_max_abs = NA_real_,
          grid_iae = NA_real_,
          grid_ise = NA_real_,
          ok = FALSE,
          error = conditionMessage(e),
          stringsAsFactors = FALSE
        ),
        eval = NULL
      )
    }
  )

  out
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

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  median(x)
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  mean(x)
}

safe_quantile <- function(x, prob) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  as.numeric(stats::quantile(x, probs = prob, names = FALSE, type = 7L))
}

summarize_runs <- function(raw_runs, raw_eval, nmulti_ref) {
  ok_runs <- raw_runs[raw_runs$ok, , drop = FALSE]
  if (nrow(ok_runs) == 0L) return(list(summary = data.frame(), drift = data.frame()))

  split_key <- interaction(ok_runs$scenario_id, ok_runs$n, ok_runs$seed, drop = TRUE)
  best_primary <- tapply(seq_len(nrow(ok_runs)), split_key, function(idx) {
    rows <- ok_runs[idx, , drop = FALSE]
    metric_name <- if (rows$family[1L] == "npreg") "pred_rmse" else "grid_iae"
    min(rows[[metric_name]], na.rm = TRUE)
  })
  ok_runs$best_primary <- as.numeric(best_primary[split_key])
  ok_runs$primary_value <- ifelse(ok_runs$family == "npreg", ok_runs$pred_rmse, ok_runs$grid_iae)
  ok_runs$within_1pct_best <- ok_runs$primary_value <= 1.01 * ok_runs$best_primary

  ref_runs <- ok_runs[ok_runs$nmulti == nmulti_ref, , drop = FALSE]
  names(ref_runs)[names(ref_runs) %in% c("elapsed_bw", "elapsed_fit", "elapsed_total", "bandwidth", "pred_rmse", "pred_mae", "pred_max_abs", "grid_iae", "grid_ise", "fval")] <-
    paste0(names(ref_runs)[names(ref_runs) %in% c("elapsed_bw", "elapsed_fit", "elapsed_total", "bandwidth", "pred_rmse", "pred_mae", "pred_max_abs", "grid_iae", "grid_ise", "fval")], "_ref")

  merged <- merge(
    ok_runs,
    ref_runs[, c("scenario_id", "n", "seed", "elapsed_bw_ref", "elapsed_fit_ref", "elapsed_total_ref", "bandwidth_ref", "pred_rmse_ref", "pred_mae_ref", "pred_max_abs_ref", "grid_iae_ref", "grid_ise_ref", "fval_ref")],
    by = c("scenario_id", "n", "seed"),
    all.x = TRUE,
    sort = FALSE
  )

  merged$bw_speedup_vs_ref <- merged$elapsed_bw_ref / merged$elapsed_bw
  merged$total_speedup_vs_ref <- merged$elapsed_total_ref / merged$elapsed_total
  merged$bandwidth_l2_vs_ref <- mapply(bandwidth_l2_diff, merged$bandwidth, merged$bandwidth_ref)
  merged$pred_rmse_delta_vs_ref <- merged$pred_rmse - merged$pred_rmse_ref
  merged$grid_iae_delta_vs_ref <- merged$grid_iae - merged$grid_iae_ref
  merged$grid_ise_delta_vs_ref <- merged$grid_ise - merged$grid_ise_ref
  merged$fval_abs_gap_vs_ref <- abs(merged$fval - merged$fval_ref)

  eval_ok <- raw_eval[raw_eval$scenario_id %in% ok_runs$scenario_id &
                        raw_eval$n %in% ok_runs$n &
                        raw_eval$seed %in% ok_runs$seed, , drop = FALSE]
  ref_eval <- eval_ok[eval_ok$nmulti == nmulti_ref, , drop = FALSE]
  names(ref_eval)[names(ref_eval) %in% c("pred", "truth")] <- c("pred_ref", "truth_ref")

  eval_merged <- merge(
    eval_ok,
    ref_eval[, c("scenario_id", "n", "seed", "eval_id", "pred_ref", "truth_ref")],
    by = c("scenario_id", "n", "seed", "eval_id"),
    all.x = TRUE,
    sort = FALSE
  )

  drift_groups <- split(eval_merged, interaction(eval_merged$scenario_id, eval_merged$n, eval_merged$seed, eval_merged$nmulti, drop = TRUE))
  drift_rows <- lapply(drift_groups, function(df) {
    data.frame(
      scenario_id = df$scenario_id[1L],
      family = df$family[1L],
      n = df$n[1L],
      seed = df$seed[1L],
      nmulti = df$nmulti[1L],
      fit_drift_rmse_vs_ref = sqrt(mean((df$pred - df$pred_ref)^2)),
      fit_drift_max_abs_vs_ref = max(abs(df$pred - df$pred_ref)),
      stringsAsFactors = FALSE
    )
  })
  drift <- do.call(rbind, drift_rows)

  merged <- merge(
    merged,
    aggregate(cbind(fit_drift_rmse_vs_ref, fit_drift_max_abs_vs_ref) ~ scenario_id + n + seed + nmulti,
              data = drift,
              FUN = safe_mean),
    by = c("scenario_id", "n", "seed", "nmulti"),
    all.x = TRUE,
    sort = FALSE
  )

  group_rows <- split(merged, interaction(merged$scenario_id, merged$n, merged$nmulti, drop = TRUE))
  summary_rows <- lapply(group_rows, function(df) {
    data.frame(
      scenario_id = df$scenario_id[1L],
      family = df$family[1L],
      n = df$n[1L],
      nmulti = df$nmulti[1L],
      runs = nrow(df),
      ok_rate = mean(df$ok),
      median_elapsed_bw = safe_median(df$elapsed_bw),
      median_elapsed_total = safe_median(df$elapsed_total),
      median_bw_speedup_vs_ref = safe_median(df$bw_speedup_vs_ref),
      median_total_speedup_vs_ref = safe_median(df$total_speedup_vs_ref),
      median_bandwidth_l2_vs_ref = safe_median(df$bandwidth_l2_vs_ref),
      p90_bandwidth_l2_vs_ref = safe_quantile(df$bandwidth_l2_vs_ref, 0.90),
      median_pred_rmse = safe_median(df$pred_rmse),
      median_pred_rmse_delta_vs_ref = safe_median(df$pred_rmse_delta_vs_ref),
      median_grid_iae = safe_median(df$grid_iae),
      median_grid_iae_delta_vs_ref = safe_median(df$grid_iae_delta_vs_ref),
      median_grid_ise_delta_vs_ref = safe_median(df$grid_ise_delta_vs_ref),
      median_fit_drift_rmse_vs_ref = safe_median(df$fit_drift_rmse_vs_ref),
      p90_fit_drift_rmse_vs_ref = safe_quantile(df$fit_drift_rmse_vs_ref, 0.90),
      median_fit_drift_max_abs_vs_ref = safe_median(df$fit_drift_max_abs_vs_ref),
      median_fval_abs_gap_vs_ref = safe_median(df$fval_abs_gap_vs_ref),
      share_within_1pct_best = mean(df$within_1pct_best),
      stringsAsFactors = FALSE
    )
  })

  list(
    summary = do.call(rbind, summary_rows),
    drift = drift
  )
}

write_report <- function(summary_df, out_file, nmulti_ref) {
  if (nrow(summary_df) == 0L) {
    writeLines(c("# nmulti default study", "", "No successful runs."), con = out_file)
    return(invisible(NULL))
  }

  split_family <- split(summary_df, summary_df$family)
  lines <- c(
    "# nmulti Default Study",
    "",
    sprintf("Reference `nmulti` for drift comparisons: `%d`.", nmulti_ref),
    "",
    "Interpretation:",
    "- `median_bw_speedup_vs_ref > 1` means the candidate is faster than the reference on bandwidth selection.",
    "- `share_within_1pct_best` is the fraction of seed-matched datasets where the candidate's truth-based primary error is within 1% of the best observed among the tested `nmulti` values.",
    "- `median_fit_drift_rmse_vs_ref` measures how much the evaluation-surface predictions move relative to the reference `nmulti`.",
    ""
  )

  for (fam in names(split_family)) {
    df <- split_family[[fam]]
    lines <- c(lines, sprintf("## %s", fam), "")
    for (sid in unique(df$scenario_id)) {
      sub <- df[df$scenario_id == sid, , drop = FALSE]
      sub <- sub[order(sub$n, sub$nmulti), , drop = FALSE]
      lines <- c(lines, sprintf("### %s", sid), "")
      hdr <- "| n | nmulti | med bw sec | med bw speedup | med fit drift RMSE | share within 1% best |"
      sep <- "|---:|---:|---:|---:|---:|---:|"
      rows <- apply(sub, 1L, function(r) {
        sprintf(
          "| %s | %s | %.3f | %.2f | %.5f | %.2f |",
          r[["n"]],
          r[["nmulti"]],
          as.numeric(r[["median_elapsed_bw"]]),
          as.numeric(r[["median_bw_speedup_vs_ref"]]),
          as.numeric(r[["median_fit_drift_rmse_vs_ref"]]),
          as.numeric(r[["share_within_1pct_best"]])
        )
      })
      lines <- c(lines, hdr, sep, rows, "")
    }
  }

  writeLines(lines, con = out_file)
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  cfg <- resolve_config(parse_args(args))
  dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
  load_backend()

  raw_runs <- list()
  raw_eval <- list()
  idx <- 0L
  total_jobs <- length(cfg$scenarios) * length(cfg$n_grid) * length(cfg$seeds) * length(cfg$nmulti_grid)

  if (cfg$show_progress) {
    cat("Running nmulti default study\n")
    cat("out_dir=", cfg$out_dir, "\n", sep = "")
    cat("scenarios=", paste(vapply(cfg$scenarios, `[[`, character(1), "id"), collapse = ","), "\n", sep = "")
    cat("n_grid=", paste(cfg$n_grid, collapse = ","), "\n", sep = "")
    cat("nmulti_grid=", paste(cfg$nmulti_grid, collapse = ","), "\n", sep = "")
    cat("seeds=", paste(cfg$seeds, collapse = ","), "\n", sep = "")
  }

  for (spec in cfg$scenarios) {
    eval_template <- build_eval_data(spec)
    for (n in cfg$n_grid) {
      for (seed in cfg$seeds) {
        for (nmulti in cfg$nmulti_grid) {
          idx <- idx + 1L
          if (cfg$show_progress) {
            cat(sprintf("[%d/%d] scenario=%s n=%d seed=%d nmulti=%d\n",
                        idx, total_jobs, spec$id, n, seed, nmulti))
          }
          res <- run_one(spec = spec, n = n, seed = seed, nmulti = nmulti, eval_template = eval_template)
          raw_runs[[length(raw_runs) + 1L]] <- res$run
          if (!is.null(res$eval)) raw_eval[[length(raw_eval) + 1L]] <- res$eval
        }
      }
    }
  }

  raw_runs_df <- do.call(rbind, raw_runs)
  raw_eval_df <- if (length(raw_eval) > 0L) do.call(rbind, raw_eval) else data.frame()
  nmulti_ref <- max(cfg$nmulti_grid)
  summaries <- summarize_runs(raw_runs_df, raw_eval_df, nmulti_ref = nmulti_ref)

  raw_runs_file <- file.path(cfg$out_dir, "nmulti_default_study_runs.csv")
  raw_eval_file <- file.path(cfg$out_dir, "nmulti_default_study_eval.csv")
  summary_file <- file.path(cfg$out_dir, "nmulti_default_study_summary.csv")
  drift_file <- file.path(cfg$out_dir, "nmulti_default_study_drift.csv")
  report_file <- file.path(cfg$out_dir, "nmulti_default_study_report.md")

  write.csv(raw_runs_df, file = raw_runs_file, row.names = FALSE)
  write.csv(raw_eval_df, file = raw_eval_file, row.names = FALSE)
  write.csv(summaries$summary, file = summary_file, row.names = FALSE)
  write.csv(summaries$drift, file = drift_file, row.names = FALSE)
  write_report(summaries$summary, out_file = report_file, nmulti_ref = nmulti_ref)

  if (cfg$show_progress) {
    cat("Runs:", raw_runs_file, "\n")
    cat("Eval:", raw_eval_file, "\n")
    cat("Summary:", summary_file, "\n")
    cat("Drift:", drift_file, "\n")
    cat("Report:", report_file, "\n")
    if (nrow(summaries$summary) > 0L) print(summaries$summary)
  }

  invisible(NULL)
}

if (sys.nframe() == 0L) main()
