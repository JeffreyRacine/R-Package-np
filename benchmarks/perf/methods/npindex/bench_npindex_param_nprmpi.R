#!/usr/bin/env Rscript

parse_args <- function(args) {
  out <- list(
    n = 100L,
    method = "ichimura",
    nmulti = 1L,
    ckertype = "gaussian",
    np_tree = FALSE,
    seed_policy = "varying",
    times = 50L,
    base_seed = 42L,
    seeds = NULL,
    nslaves = 1L,
    out_raw = "/tmp/npindex_mpi_bench_raw.csv",
    out_summary = "/tmp/npindex_mpi_bench_summary.csv",
    show_progress = TRUE
  )

  if (length(args) == 0L) return(out)
  for (a in args) {
    if (!grepl("^--", a)) stop("Bad arg: ", a)
    raw <- sub("^--", "", a)
    key <- sub("=.*$", "", raw)
    val <- if (grepl("=", raw, fixed = TRUE)) sub("^[^=]*=", "", raw) else "TRUE"

    if (key == "n") out$n <- as.integer(val)
    else if (key == "method") out$method <- val
    else if (key == "nmulti") out$nmulti <- as.integer(val)
    else if (key == "ckertype" || key == "ckerpair") out$ckertype <- val
    else if (key == "np_tree" || key == "np.trees" || key == "np_trees") out$np_tree <- as.logical(val)
    else if (key == "seed_policy") out$seed_policy <- val
    else if (key == "times") out$times <- as.integer(val)
    else if (key == "base_seed") out$base_seed <- as.integer(val)
    else if (key == "seeds") out$seeds <- as.integer(strsplit(val, ",", fixed = TRUE)[[1]])
    else if (key == "nslaves" || key == "rslaves") out$nslaves <- as.integer(val)
    else if (key == "out_raw") out$out_raw <- val
    else if (key == "out_summary") out$out_summary <- val
    else if (key == "show_progress") out$show_progress <- as.logical(val)
    else stop("Unknown arg: ", key)
  }

  if (!out$method %in% c("ichimura", "kleinspady")) stop("method must be ichimura or kleinspady")
  if (!is.finite(out$nmulti) || out$nmulti < 1L) stop("nmulti must be an integer >= 1")
  if (!out$ckertype %in% c("gaussian", "epanechnikov")) stop("ckertype must be gaussian or epanechnikov")
  if (!out$seed_policy %in% c("fixed", "varying")) stop("seed_policy must be fixed or varying")
  if (!is.finite(out$nslaves) || out$nslaves < 1L) stop("nslaves must be >= 1")

  out
}

make_seeds <- function(cfg) {
  if (!is.null(cfg$seeds)) return(cfg$seeds)
  cfg$base_seed
}

run_case <- function(seed, cfg) {
  first_num_or_na <- function(x) {
    if (is.null(x) || length(x) == 0L) return(NA_real_)
    as.numeric(x)[1]
  }
  count_or_na <- function(x) {
    if (is.null(x)) return(NA_integer_)
    nr <- nrow(x)
    if (!is.null(nr) && length(nr) > 0L) return(as.integer(nr)[1])
    as.integer(length(x))[1]
  }
  bw_str <- function(bw) {
    vals <- NULL
    if (!is.null(bw$bw) && length(bw$bw) > 0L) vals <- suppressWarnings(as.numeric(unlist(bw$bw)))
    if ((is.null(vals) || length(vals) == 0L) && !is.null(bw$bandwidth) && length(bw$bandwidth) > 0L) vals <- suppressWarnings(as.numeric(unlist(bw$bandwidth)))
    vals <- vals[is.finite(vals)]
    if (is.null(vals) || length(vals) == 0L) return("")
    paste(signif(vals, 8), collapse = ";")
  }
  first_num_multi <- function(...) {
    xs <- list(...)
    for (x in xs) if (!is.null(x) && length(x) > 0L) return(as.numeric(x)[1])
    NA_real_
  }

  iter_seeds <- if (cfg$seed_policy == "fixed") rep(seed, cfg$times) else seed + seq_len(cfg$times) - 1L
  res <- vector("list", cfg$times)
  idx <- 0L

  mb <- microbenchmark::microbenchmark(
    run = {
      idx <- idx + 1L
      s <- iter_seeds[[idx]]
      out <- tryCatch({
        t0 <- proc.time()[["elapsed"]]
        eval(substitute(
          {
            set.seed(S)
            x1 <- runif(N, min = -1, max = 1)
            x2 <- runif(N, min = -1, max = 1)
            if (MTH == "ichimura") {
              y <- x1 - x2 + rnorm(N)
            } else {
              y <- ifelse(x1 + x2 + rnorm(N) > 0, 1L, 0L)
            }
            dat <- data.frame(y = y, x1 = x1, x2 = x2)
            bw <- npindexbw(y ~ x1 + x2,
                            method = MTH,
                            nmulti = NM,
                            ckertype = CK,
                            data = dat)
          },
          list(S = s, N = cfg$n, MTH = cfg$method, NM = cfg$nmulti, CK = cfg$ckertype)
        ))
        bw_elapsed <- proc.time()[["elapsed"]] - t0

        t1 <- proc.time()[["elapsed"]]
        eval(substitute(
          {
            fit <- npindex(bws = bw, data = dat)
          },
          list()
        ))
        fit_elapsed <- proc.time()[["elapsed"]] - t1

        fit_metric <- if (cfg$method == "ichimura") first_num_or_na(fit$MSE) else first_num_or_na(fit$fit.mcfadden)

        list(
          ok = !is.null(fit),
          elapsed_bw = bw_elapsed,
          elapsed_fit = fit_elapsed,
          fval = first_num_or_na(bw$fval),
          ifval = first_num_or_na(bw$ifval),
          eval_count = count_or_na(bw$eval.history),
          invalid_count = count_or_na(bw$invalid.history),
          num_fval = first_num_multi(bw$num.fval, bw$num.feval),
          log_likelihood = fit_metric,
          bw = bw_str(bw),
          error = ""
        )
      }, error = function(e) {
        list(
          ok = FALSE,
          elapsed_bw = NA_real_,
          elapsed_fit = NA_real_,
          fval = NA_real_,
          ifval = NA_real_,
          eval_count = NA_integer_,
          invalid_count = NA_integer_,
          num_fval = NA_real_,
          log_likelihood = NA_real_,
          bw = "",
          error = conditionMessage(e)
        )
      })
      res[[idx]] <- out
      invisible(NULL)
    },
    times = cfg$times
  )

  total_secs <- as.numeric(mb$time) / 1e9

  data.frame(
    backend = rep("npRmpi", cfg$times),
    n = rep(cfg$n, cfg$times),
    bwmethod = rep(cfg$method, cfg$times),
    method = rep(cfg$method, cfg$times),
    nmulti = rep(cfg$nmulti, cfg$times),
    ckertype = rep(cfg$ckertype, cfg$times),
    np_tree = rep(cfg$np_tree, cfg$times),
    seed_policy = rep(cfg$seed_policy, cfg$times),
    seed = iter_seeds,
    iter = seq_len(cfg$times),
    ok = vapply(res, function(x) isTRUE(x$ok), logical(1)),
    elapsed_bw = vapply(res, function(x) x$elapsed_bw, numeric(1)),
    elapsed_fit = vapply(res, function(x) x$elapsed_fit, numeric(1)),
    elapsed_total = total_secs,
    fval = vapply(res, function(x) x$fval, numeric(1)),
    ifval = vapply(res, function(x) x$ifval, numeric(1)),
    eval_count = vapply(res, function(x) x$eval_count, integer(1)),
    invalid_count = vapply(res, function(x) x$invalid_count, integer(1)),
    num_fval = vapply(res, function(x) x$num_fval, numeric(1)),
    log_likelihood = vapply(res, function(x) x$log_likelihood, numeric(1)),
    bw = vapply(res, function(x) x$bw, character(1)),
    error = vapply(res, function(x) x$error, character(1)),
    stringsAsFactors = FALSE
  )
}

summarize_results <- function(df) {
  okdf <- df[df$ok, , drop = FALSE]
  if (nrow(okdf) == 0L) return(data.frame())

  data.frame(
    backend = okdf$backend[1],
    n = okdf$n[1],
    bwmethod = okdf$bwmethod[1],
    method = okdf$method[1],
    nmulti = okdf$nmulti[1],
    ckertype = okdf$ckertype[1],
    np_tree = okdf$np_tree[1],
    seed_policy = okdf$seed_policy[1],
    runs = nrow(okdf),
    mean_elapsed_total = mean(okdf$elapsed_total),
    median_elapsed_total = median(okdf$elapsed_total),
    mean_elapsed_bw = mean(okdf$elapsed_bw),
    median_elapsed_bw = median(okdf$elapsed_bw),
    mean_elapsed_fit = mean(okdf$elapsed_fit),
    median_elapsed_fit = median(okdf$elapsed_fit),
    mean_fval = mean(okdf$fval),
    median_fval = median(okdf$fval),
    mean_eval_count = mean(okdf$eval_count),
    median_eval_count = median(okdf$eval_count),
    mean_invalid_count = mean(okdf$invalid_count),
    median_invalid_count = median(okdf$invalid_count),
    mean_num_fval = mean(okdf$num_fval),
    median_num_fval = median(okdf$num_fval),
    mean_log_likelihood = mean(okdf$log_likelihood),
    median_log_likelihood = median(okdf$log_likelihood),
    stringsAsFactors = FALSE
  )
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  cfg <- parse_args(args)

  suppressPackageStartupMessages(library(npRmpi))
  suppressPackageStartupMessages(library(microbenchmark))

  npRmpi.start(nslaves = cfg$nslaves)
  on.exit(try(npRmpi.stop(force = TRUE), silent = TRUE), add = TRUE)

  options(npRmpi.autodispatch = TRUE, np.messages = FALSE, np.tree = cfg)
  eval(substitute(
    options(npRmpi.autodispatch = TRUE, np.messages = FALSE, np.tree = TF),
    list(TF = cfg$np_tree)
  ))

  seeds <- make_seeds(cfg)

  if (cfg$show_progress) {
    cat("Running npindex npRmpi benchmark with", length(seeds), "seeds\n")
    cat("backend=npRmpi",
        " n=", cfg$n,
        " nslaves=", cfg$nslaves,
        " method=", cfg$method,
        " nmulti=", cfg$nmulti,
        " ckertype=", cfg$ckertype,
        " np_tree=", cfg$np_tree,
        " seed_policy=", cfg$seed_policy, "\n", sep = "")
  }

  rows <- vector("list", length(seeds))
  for (i in seq_along(seeds)) {
    if (cfg$show_progress) cat("  run", i, "seed", seeds[[i]], "\n")
    rows[[i]] <- run_case(seeds[[i]], cfg)
  }

  raw <- do.call(rbind, rows)
  summary <- summarize_results(raw)

  write.table(raw,
              file = cfg$out_raw,
              sep = ",",
              row.names = FALSE,
              col.names = !file.exists(cfg$out_raw),
              append = file.exists(cfg$out_raw),
              qmethod = "double")

  if (nrow(summary) > 0L) {
    write.table(summary,
                file = cfg$out_summary,
                sep = ",",
                row.names = FALSE,
                col.names = !file.exists(cfg$out_summary),
                append = file.exists(cfg$out_summary),
                qmethod = "double")
  }

  if (cfg$show_progress) {
    cat("Raw:", cfg$out_raw, "\n")
    cat("Summary:", cfg$out_summary, "\n")
    if (nrow(summary) > 0L) print(summary)
  }
}

if (sys.nframe() == 0) main()
