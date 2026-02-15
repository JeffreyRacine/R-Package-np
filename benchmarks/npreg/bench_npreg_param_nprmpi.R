#!/usr/bin/env Rscript

parse_args <- function(args) {
  out <- list(
    n = 100L,
    regtype = "lc",
    bwmethod = "cv.ls",
    nmulti = 1L,
    ckertype = "gaussian",
    np_tree = FALSE,
    seed_policy = "varying",
    times = 5L,
    base_seed = 42L,
    seeds = NULL,
    nslaves = 1L,
    out_raw = "/tmp/npreg_nprmpi_bench_raw.csv",
    out_summary = "/tmp/npreg_nprmpi_bench_summary.csv",
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
    else if (key == "bwmethod") out$bwmethod <- val
    else if (key == "nmulti") out$nmulti <- as.integer(val)
    else if (key == "ckertype") out$ckertype <- val
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

  if (!out$regtype %in% c("lc", "ll")) stop("regtype must be lc or ll")
  if (!out$bwmethod %in% c("cv.ls", "cv.aic")) stop("bwmethod must be cv.ls or cv.aic")
  if (!is.finite(out$nmulti) || out$nmulti < 1L) stop("nmulti must be an integer >= 1")
  if (!out$ckertype %in% c("gaussian", "epanechnikov")) stop("ckertype must be gaussian or epanechnikov")
  if (!out$seed_policy %in% c("fixed", "varying")) stop("seed_policy must be fixed or varying")

  out
}

make_seeds <- function(cfg) {
  if (!is.null(cfg$seeds)) return(cfg$seeds)
  if (cfg$seed_policy == "fixed") rep(cfg$base_seed, cfg$times)
  else cfg$base_seed + seq_len(cfg$times) - 1L
}

mk_data <- function(seed, n) {
  set.seed(seed)
  x1 <- runif(n)
  x2 <- runif(n)
  z1 <- factor(rbinom(n, 3, .5))
  z2 <- ordered(rbinom(n, 3, .5))
  y <- cos(2 * pi * x1) + as.numeric(as.character(z1)) + rnorm(n, sd = .25)
  data.frame(y, x1, x2, z1, z2)
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

  d <- mk_data(seed, cfg$n)
  mpi.bcast.Robj2slave(d)

  err <- ""
  ok <- FALSE

  bw <- NULL
  fit <- NULL
  elapsed_bw <- NA_real_
  elapsed_fit <- NA_real_
  elapsed_total <- NA_real_
  num_fval <- NA_real_

  t0 <- proc.time()[["elapsed"]]
  t_bw <- system.time({
    tryCatch({
      bw <- eval(substitute(
        mpi.bcast.cmd(
          bw <- npregbw(y ~ x1 + x2 + z1 + z2,
                        regtype = REGTYPE,
                        bwmethod = BWMETHOD,
                        nmulti = NMULTI,
                        ckertype = CKERTYPE,
                        data = d),
          caller.execute = TRUE
        ),
        list(REGTYPE = cfg$regtype, BWMETHOD = cfg$bwmethod, NMULTI = cfg$nmulti, CKERTYPE = cfg$ckertype)
      ))
    }, error = function(e) err <<- conditionMessage(e))
  })
  elapsed_bw <- unname(t_bw[["elapsed"]])

  t_fit <- system.time({
    if (!is.null(bw)) {
      tryCatch({
        fit <- mpi.bcast.cmd(npreg(bws = bw, data = d), caller.execute = TRUE)
      }, error = function(e) err <<- conditionMessage(e))
    }
  })
  elapsed_fit <- unname(t_fit[["elapsed"]])
  elapsed_total <- proc.time()[["elapsed"]] - t0
  ok <- !is.null(bw) && !is.null(fit) && identical(err, "")
  if (!is.null(bw)) num_fval <- first_num_or_na(bw$num.fval)

  data.frame(
    backend = "npRmpi",
    n = cfg$n,
    regtype = cfg$regtype,
    bwmethod = cfg$bwmethod,
    nmulti = cfg$nmulti,
    ckertype = cfg$ckertype,
    np_tree = cfg$np_tree,
    seed_policy = cfg$seed_policy,
    seed = seed,
    ok = ok,
    elapsed_bw = elapsed_bw,
    elapsed_fit = elapsed_fit,
    elapsed_total = elapsed_total,
    fval = if (!is.null(bw)) first_num_or_na(bw$fval) else NA_real_,
    ifval = if (!is.null(bw)) first_num_or_na(bw$ifval) else NA_real_,
    eval_count = if (!is.null(bw)) count_or_na(bw$eval.history) else NA_integer_,
    invalid_count = if (!is.null(bw)) count_or_na(bw$invalid.history) else NA_integer_,
    num_fval = num_fval,
    bw = if (!is.null(bw)) bw_str_or_empty(bw$bw) else "",
    error = err,
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

  suppressPackageStartupMessages(library(npRmpi))
  npRmpi.start(nslaves = cfg$nslaves)
  on.exit(try(npRmpi.stop(force = TRUE), silent = TRUE), add = TRUE)

  options(np.messages = FALSE, np.tree = cfg$np_tree)
  eval(substitute(
    mpi.bcast.cmd(options(np.messages = FALSE, np.tree = TF), caller.execute = TRUE),
    list(TF = cfg$np_tree)
  ))

  seeds <- make_seeds(cfg)

  if (cfg$show_progress) {
    cat("Running npreg npRmpi benchmark with", length(seeds), "seeds\n")
    cat("backend=npRmpi",
        " n=", cfg$n,
        " nslaves=", cfg$nslaves,
        " regtype=", cfg$regtype,
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
