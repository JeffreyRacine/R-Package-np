#!/usr/bin/env Rscript

parse_args <- function(args) {
  parse_degree <- function(val) {
    parts <- strsplit(val, ",", fixed = TRUE)[[1]]
    out <- as.integer(parts)
    if (length(out) == 0L || anyNA(out)) stop("degree must be comma-separated non-negative integers (e.g. 2,2)")
    if (any(out < 0L)) stop("degree must contain non-negative integers")
    out
  }

  out <- list(
    n = 100L,
    regtype = "ll",
    basis = "glp",
    degree = c(2L, 2L),
    bernstein.basis = FALSE,
    bwmethod = "cv.ls",
    nmulti = 1L,
    ckertype = "gaussian",
    np_tree = FALSE,
    seed_policy = "varying",
    times = 50L,
    base_seed = 42L,
    seeds = NULL,
    out_raw = "/tmp/npreg_serial_bench_raw.csv",
    out_summary = "/tmp/npreg_serial_bench_summary.csv",
    show_progress = TRUE
  )

  if (length(args) == 0L) return(out)
  for (a in args) {
    if (!grepl("^--", a)) stop("Bad arg: ", a)
    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
    key <- kv[1]
    val <- if (length(kv) > 1L) kv[2] else "TRUE"

    if (key == "n") out$n <- as.integer(val)
    else if (key == "regtype") out$regtype <- val
    else if (key == "basis") out$basis <- val
    else if (key == "degree") out$degree <- parse_degree(val)
    else if (key == "bernstein.basis" || key == "bernstein_basis") out$bernstein.basis <- as.logical(val)
    else if (key == "bwmethod") out$bwmethod <- val
    else if (key == "nmulti") out$nmulti <- as.integer(val)
    else if (key == "ckertype") out$ckertype <- val
    else if (key == "np_tree" || key == "np.trees" || key == "np_trees") out$np_tree <- as.logical(val)
    else if (key == "seed_policy") out$seed_policy <- val
    else if (key == "times") out$times <- as.integer(val)
    else if (key == "base_seed") out$base_seed <- as.integer(val)
    else if (key == "seeds") out$seeds <- as.integer(strsplit(val, ",", fixed = TRUE)[[1]])
    else if (key == "out_raw") out$out_raw <- val
    else if (key == "out_summary") out$out_summary <- val
    else if (key == "show_progress") out$show_progress <- as.logical(val)
    else stop("Unknown arg: ", key)
  }

  if (!out$regtype %in% c("lc", "ll", "lp")) stop("regtype must be lc, ll, or lp")
  if (!out$basis %in% c("glp", "additive", "tensor")) stop("basis must be glp, additive, or tensor")
  if (!out$bwmethod %in% c("cv.ls", "cv.aic")) stop("bwmethod must be cv.ls or cv.aic")
  if (!is.finite(out$nmulti) || out$nmulti < 1L) stop("nmulti must be an integer >= 1")
  if (!out$ckertype %in% c("gaussian", "epanechnikov")) stop("ckertype must be gaussian or epanechnikov")
  if (!out$seed_policy %in% c("fixed", "varying")) stop("seed_policy must be fixed or varying")
  if (out$regtype == "lp" && length(out$degree) < 1L) stop("degree must be supplied for regtype=lp")

  out
}

make_seeds <- function(cfg) {
  if (!is.null(cfg$seeds)) return(cfg$seeds)
  cfg$base_seed
}

load_backend <- function() {
  suppressPackageStartupMessages(library(np))
  suppressPackageStartupMessages(library(microbenchmark))
  options(np.messages = FALSE)
}

mk_data <- function(seed, n) {
  set.seed(seed)
  x1 <- runif(n)
  x2 <- runif(n)
  z1 <- factor(rbinom(n, 3, .5))
  z2 <- ordered(rbinom(n, 3, .5))
  y <- cos(2 * pi * x1) + as.numeric(as.character(z1)) + rnorm(n, sd = .25)
  data.frame(y = y, x1 = x1, x2 = x2, z1 = z1, z2 = z2)
}

run_case <- function(seed, cfg) {
  first_num_or_na <- function(x) {
    if (is.null(x) || length(x) == 0L) return(NA_real_)
    as.numeric(x)[1]
  }
  count_or_na <- function(x) {
    if (is.null(x)) return(NA_integer_)
    nr <- nrow(x)
    if (is.null(nr) || length(nr) == 0L) return(NA_integer_)
    as.integer(nr)[1]
  }
  bw_str_or_empty <- function(x) {
    if (is.null(x) || length(x) == 0L) return("")
    paste(signif(as.numeric(x), 8), collapse = ";")
  }
  degree_str <- paste(cfg$degree, collapse = ",")

  iter_seeds <- if (cfg$seed_policy == "fixed") rep(seed, cfg$times) else seed + seq_len(cfg$times) - 1L
  res <- vector("list", cfg$times)
  idx <- 0L

  mb <- microbenchmark::microbenchmark(
    run = {
      idx <- idx + 1L
      s <- iter_seeds[[idx]]
      out <- tryCatch({
        dat <- mk_data(s, cfg$n)

        t0 <- proc.time()[["elapsed"]]
        bw_args <- list(
          formula = y ~ x1 + x2 + z1 + z2,
          regtype = cfg$regtype,
          bwmethod = cfg$bwmethod,
          nmulti = cfg$nmulti,
          ckertype = cfg$ckertype,
          data = dat
        )
        if (cfg$regtype == "lp") {
          bw_args$basis <- cfg$basis
          bw_args$degree <- cfg$degree
          bw_args$bernstein.basis <- isTRUE(cfg$bernstein.basis)
        }
        bw <- do.call(npregbw, bw_args)
        bw_elapsed <- proc.time()[["elapsed"]] - t0

        t1 <- proc.time()[["elapsed"]]
        fit <- npreg(bws = bw, data = dat)
        fit_elapsed <- proc.time()[["elapsed"]] - t1

        list(
          ok = !is.null(fit),
          elapsed_bw = bw_elapsed,
          elapsed_fit = fit_elapsed,
          fval = first_num_or_na(bw$fval),
          ifval = first_num_or_na(bw$ifval),
          eval_count = count_or_na(bw$eval.history),
          invalid_count = count_or_na(bw$invalid.history),
          num_fval = first_num_or_na(bw$num.fval),
          bw = bw_str_or_empty(bw$bw),
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
    backend = rep("np", cfg$times),
    n = rep(cfg$n, cfg$times),
    regtype = rep(cfg$regtype, cfg$times),
    basis = rep(if (cfg$regtype == "lp") cfg$basis else NA_character_, cfg$times),
    degree = rep(if (cfg$regtype == "lp") degree_str else "", cfg$times),
    bernstein.basis = rep(if (cfg$regtype == "lp") isTRUE(cfg$bernstein.basis) else NA, cfg$times),
    bwmethod = rep(cfg$bwmethod, cfg$times),
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
    regtype = okdf$regtype[1],
    basis = if ("basis" %in% names(okdf)) okdf$basis[1] else NA_character_,
    degree = if ("degree" %in% names(okdf)) okdf$degree[1] else "",
    bernstein.basis = if ("bernstein.basis" %in% names(okdf)) okdf$bernstein.basis[1] else NA,
    bwmethod = okdf$bwmethod[1],
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
    stringsAsFactors = FALSE
  )
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  cfg <- parse_args(args)
  load_backend()
  options(np.tree = cfg$np_tree)
  seeds <- make_seeds(cfg)

  if (cfg$show_progress) {
    cat("Running npreg serial benchmark with", length(seeds), "seeds\n")
    cat("backend=np",
        " n=", cfg$n,
        " regtype=", cfg$regtype,
        if (cfg$regtype == "lp") paste0(" basis=", cfg$basis) else "",
        if (cfg$regtype == "lp") paste0(" degree=", paste(cfg$degree, collapse = ",")) else "",
        if (cfg$regtype == "lp") paste0(" bernstein.basis=", isTRUE(cfg$bernstein.basis)) else "",
        " bwmethod=", cfg$bwmethod,
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
