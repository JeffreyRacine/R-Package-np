#!/usr/bin/env Rscript

parse_args <- function(args) {
  out <- list(
    fun = "",
    n = 1000L,
    nslaves = 1L,
    times = 50L,
    seed_policy = "varying",
    base_seed = 42L,
    seeds = NULL,
    out_raw = "/tmp/nprmpi_oneoff_bench_raw.csv",
    out_summary = "/tmp/nprmpi_oneoff_bench_summary.csv",
    show_progress = TRUE
  )

  if (length(args) == 0L) return(out)
  for (a in args) {
    if (!grepl("^--", a)) stop("Bad arg: ", a)
    raw <- sub("^--", "", a)
    key <- sub("=.*$", "", raw)
    val <- if (grepl("=", raw, fixed = TRUE)) sub("^[^=]*=", "", raw) else "TRUE"

    if (key == "fun") out$fun <- val
    else if (key == "n") out$n <- as.integer(val)
    else if (key == "nslaves" || key == "rslaves") out$nslaves <- as.integer(val)
    else if (key == "times") out$times <- as.integer(val)
    else if (key == "seed_policy") out$seed_policy <- val
    else if (key == "base_seed") out$base_seed <- as.integer(val)
    else if (key == "seeds") out$seeds <- as.integer(strsplit(val, ",", fixed = TRUE)[[1]])
    else if (key == "out_raw") out$out_raw <- val
    else if (key == "out_summary") out$out_summary <- val
    else if (key == "show_progress") out$show_progress <- as.logical(val)
    else stop("Unknown arg: ", key)
  }

  valid <- c("npcmstest", "npconmode", "npcopula", "npdeneqtest", "npdeptest", "npqreg",
             "npregiv", "npsdeptest", "npsigtest", "npsymtest", "npunitest")
  if (!out$fun %in% valid) stop("--fun must be one of: ", paste(valid, collapse = ", "))
  if (!out$seed_policy %in% c("fixed", "varying")) stop("seed_policy must be fixed or varying")
  if (!is.finite(out$n) || out$n < 10L) stop("n must be >= 10")
  if (!is.finite(out$times) || out$times < 1L) stop("times must be >= 1")
  if (!is.finite(out$nslaves) || out$nslaves < 1L) stop("nslaves must be >= 1")

  out
}

is_real_data_fun <- function(fun) {
  fun %in% c("npcmstest", "npconmode", "npqreg")
}

make_seeds <- function(cfg) {
  if (!is.null(cfg$seeds)) return(cfg$seeds)
  cfg$base_seed
}

sig_numeric <- function(x, k = 25L) {
  nums <- suppressWarnings(as.numeric(unlist(x)))
  nums <- nums[is.finite(nums)]
  if (length(nums) == 0L) return(NA_real_)
  mean(nums[seq_len(min(k, length(nums)))])
}

first_num_or_na <- function(x) {
  if (is.null(x) || length(x) == 0L) return(NA_real_)
  suppressWarnings(as.numeric(x)[1])
}

run_one <- function(fun, n, seed) {
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)

  out <- list(
    ok = FALSE,
    n_actual = NA_integer_,
    elapsed_setup = NA_real_,
    elapsed_main = NA_real_,
    fval = NA_real_,
    sig_main = NA_real_,
    sig_aux = NA_real_,
    error = ""
  )

  tryCatch({
    if (fun == "npcmstest") {
      data("oecdpanel", package = "npRmpi")
      d <- oecdpanel
      d$oecd <- factor(d$oecd)
      d$year <- factor(d$year)
      out$n_actual <- nrow(d)

      model <- lm(growth ~ oecd + year + initgdp + I(initgdp^2) + I(initgdp^3) +
                    I(initgdp^4) + popgro + inv + humancap + I(humancap^2) + I(humancap^3) - 1,
                  data = d, x = TRUE, y = TRUE)
      X <- data.frame(d$oecd, d$year, d$initgdp, d$popgro, d$inv, d$humancap)
      names(X) <- c("oecd", "year", "initgdp", "popgro", "inv", "humancap")
      ydat <- d$growth

      t1 <- proc.time()[["elapsed"]]
      ans <- npcmstest(model = model, xdat = X, ydat = ydat, nmulti = 1L, ftol = 0.01, tol = 0.01)
      out$elapsed_main <- proc.time()[["elapsed"]] - t1
      out$sig_main <- sig_numeric(ans)
    } else if (fun == "npconmode") {
      library(MASS)
      data("birthwt", package = "MASS")
      d <- birthwt
      d$low <- factor(d$low)
      d$smoke <- factor(d$smoke)
      d$race <- factor(d$race)
      d$ht <- factor(d$ht)
      d$ui <- factor(d$ui)
      d$ftv <- ordered(d$ftv)
      out$n_actual <- nrow(d)

      t0 <- proc.time()[["elapsed"]]
      bw <- npcdensbw(low ~ smoke + race + ht + ui + ftv + age + lwt, data = d)
      out$elapsed_setup <- proc.time()[["elapsed"]] - t0

      t1 <- proc.time()[["elapsed"]]
      ans <- npconmode(bws = bw)
      out$elapsed_main <- proc.time()[["elapsed"]] - t1
      out$fval <- first_num_or_na(bw$fval)
      out$sig_aux <- sig_numeric(bw$bw)
      out$sig_main <- sig_numeric(ans)
    } else if (fun == "npcopula") {
      set.seed(seed)
      library(MASS)
      out$n_actual <- as.integer(n)
      rho <- 0.95
      mu <- c(0, 0)
      Sigma <- matrix(c(1, rho, rho, 1), 2, 2)
      mydat <- mvrnorm(n = n, mu = mu, Sigma = Sigma)
      mydat <- data.frame(x = mydat[, 1], y = mydat[, 2])
      grid.seq <- seq(0, 1, length.out = 25L)
      grid.dat <- cbind(grid.seq, grid.seq)

      t0 <- proc.time()[["elapsed"]]
      bw <- npudistbw(~ x + y, data = mydat)
      out$elapsed_setup <- proc.time()[["elapsed"]] - t0

      t1 <- proc.time()[["elapsed"]]
      ans <- npcopula(bws = bw, data = mydat, u = grid.dat)
      out$elapsed_main <- proc.time()[["elapsed"]] - t1
      out$fval <- first_num_or_na(bw$fval)
      out$sig_aux <- sig_numeric(bw$bw)
      out$sig_main <- sig_numeric(ans$copula)
    } else if (fun == "npdeneqtest") {
      set.seed(seed)
      out$n_actual <- as.integer(n)
      sample.A <- data.frame(x = rnorm(n))
      sample.B <- data.frame(x = rnorm(n))
      t1 <- proc.time()[["elapsed"]]
      ans <- npdeneqtest(sample.A, sample.B, boot.num = 99)
      out$elapsed_main <- proc.time()[["elapsed"]] - t1
      out$sig_main <- sig_numeric(ans)
    } else if (fun == "npdeptest") {
      set.seed(seed)
      out$n_actual <- as.integer(n)
      x <- rnorm(n)
      y <- 1 + x + rnorm(n)
      y.fit <- fitted(lm(y ~ x))
      t1 <- proc.time()[["elapsed"]]
      ans <- npdeptest(y, y.fit, boot.num = 99, method = "summation")
      out$elapsed_main <- proc.time()[["elapsed"]] - t1
      out$sig_main <- sig_numeric(ans)
    } else if (fun == "npqreg") {
      data("Italy", package = "npRmpi")
      d <- Italy
      out$n_actual <- nrow(d)
      t0 <- proc.time()[["elapsed"]]
      bw <- npcdistbw(gdp ~ ordered(year), data = d)
      out$elapsed_setup <- proc.time()[["elapsed"]] - t0

      t1 <- proc.time()[["elapsed"]]
      m25 <- npqreg(bws = bw, tau = 0.25)
      m50 <- npqreg(bws = bw, tau = 0.50)
      m75 <- npqreg(bws = bw, tau = 0.75)
      out$elapsed_main <- proc.time()[["elapsed"]] - t1
      out$fval <- first_num_or_na(bw$fval)
      out$sig_aux <- sig_numeric(bw$bw)
      out$sig_main <- sig_numeric(list(m25, m50, m75))
    } else if (fun == "npregiv") {
      set.seed(seed)
      out$n_actual <- as.integer(n)
      v <- rnorm(n, mean = 0, sd = 0.27)
      eps <- rnorm(n, mean = 0, sd = 0.05)
      u <- -0.5 * v + eps
      w <- rnorm(n, mean = 0, sd = 1)
      z <- 0.2 * w + v
      y <- z^2 + u
      t1 <- proc.time()[["elapsed"]]
      ans <- npregiv(y = y, z = z, w = w)
      out$elapsed_main <- proc.time()[["elapsed"]] - t1
      out$sig_main <- sig_numeric(ans)
    } else if (fun == "npsdeptest") {
      set.seed(seed)
      out$n_actual <- as.integer(n)
      ar.series <- function(phi, epsilon) {
        m <- length(epsilon)
        series <- numeric(m)
        series[1] <- epsilon[1] / (1 - phi)
        for (i in 2:m) series[i] <- phi * series[i - 1] + epsilon[i]
        series
      }
      yt <- ar.series(0.95, rnorm(n))
      t1 <- proc.time()[["elapsed"]]
      ans <- npsdeptest(yt, lag.num = 2, boot.num = 399, method = "summation")
      out$elapsed_main <- proc.time()[["elapsed"]] - t1
      out$sig_main <- sig_numeric(ans)
    } else if (fun == "npsigtest") {
      options(npRmpi.autodispatch = TRUE)
      set.seed(seed)
      out$n_actual <- as.integer(n)
      z <- factor(rbinom(n, 1, 0.5))
      x1 <- rnorm(n)
      x2 <- runif(n, -2, 2)
      y <- x1 + x2 + rnorm(n)
      d <- data.frame(z = z, x1 = x1, x2 = x2, y = y)

      t0 <- proc.time()[["elapsed"]]
      model <- npreg(y ~ z + x1 + x2, regtype = "ll", bwmethod = "cv.aic", data = d)
      out$elapsed_setup <- proc.time()[["elapsed"]] - t0

      t1 <- proc.time()[["elapsed"]]
      ans <- npsigtest(model)
      out$elapsed_main <- proc.time()[["elapsed"]] - t1
      out$sig_main <- sig_numeric(ans)
      out$sig_aux <- sig_numeric(model)
    } else if (fun == "npsymtest") {
      set.seed(seed)
      out$n_actual <- as.integer(n)
      ar.series <- function(phi, epsilon) {
        m <- length(epsilon)
        series <- numeric(m)
        series[1] <- epsilon[1] / (1 - phi)
        for (i in 2:m) series[i] <- phi * series[i - 1] + epsilon[i]
        series
      }
      yt <- ar.series(0.5, rnorm(n))
      t1 <- proc.time()[["elapsed"]]
      ans <- npsymtest(yt, boot.num = 399, boot.method = "geom", method = "summation")
      out$elapsed_main <- proc.time()[["elapsed"]] - t1
      out$sig_main <- sig_numeric(ans)
    } else if (fun == "npunitest") {
      set.seed(seed)
      out$n_actual <- as.integer(n)
      x <- rnorm(n)
      y <- rnorm(n)
      t1 <- proc.time()[["elapsed"]]
      ans <- npunitest(x, y, method = "summation", bootstrap = TRUE)
      out$elapsed_main <- proc.time()[["elapsed"]] - t1
      out$sig_main <- sig_numeric(ans)
    }

    out$ok <- TRUE
    out
  }, error = function(e) {
    out$error <- conditionMessage(e)
    out
  })
}

run_seed <- function(cfg, seed) {
  iter_seeds <- if (cfg$seed_policy == "fixed") rep(seed, cfg$times) else seed + seq_len(cfg$times) - 1L
  res <- vector("list", cfg$times)
  idx <- 0L

  mb <- microbenchmark::microbenchmark(
    run = {
      idx <- idx + 1L
      res[[idx]] <- run_one(cfg$fun, cfg$n, iter_seeds[[idx]])
      invisible(NULL)
    },
    times = cfg$times
  )

  data.frame(
    backend = rep("npRmpi", cfg$times),
    np_fun = rep(cfg$fun, cfg$times),
    n_requested = rep(cfg$n, cfg$times),
    n_actual = vapply(res, function(x) x$n_actual, integer(1)),
    nslaves = rep(cfg$nslaves, cfg$times),
    seed_policy = rep(cfg$seed_policy, cfg$times),
    seed = iter_seeds,
    iter = seq_len(cfg$times),
    ok = vapply(res, function(x) isTRUE(x$ok), logical(1)),
    elapsed_setup = vapply(res, function(x) x$elapsed_setup, numeric(1)),
    elapsed_main = vapply(res, function(x) x$elapsed_main, numeric(1)),
    elapsed_total = as.numeric(mb$time) / 1e9,
    fval = vapply(res, function(x) x$fval, numeric(1)),
    sig_main = vapply(res, function(x) x$sig_main, numeric(1)),
    sig_aux = vapply(res, function(x) x$sig_aux, numeric(1)),
    error = vapply(res, function(x) x$error, character(1)),
    stringsAsFactors = FALSE
  )
}

summarize_results <- function(df) {
  okdf <- df[df$ok, , drop = FALSE]
  if (nrow(okdf) == 0L) return(data.frame())

  data.frame(
    backend = okdf$backend[1],
    np_fun = okdf$np_fun[1],
    n_requested = okdf$n_requested[1],
    n_actual = okdf$n_actual[1],
    nslaves = okdf$nslaves[1],
    seed_policy = okdf$seed_policy[1],
    runs = nrow(okdf),
    mean_elapsed_setup = mean(okdf$elapsed_setup, na.rm = TRUE),
    median_elapsed_setup = median(okdf$elapsed_setup, na.rm = TRUE),
    mean_elapsed_main = mean(okdf$elapsed_main, na.rm = TRUE),
    median_elapsed_main = median(okdf$elapsed_main, na.rm = TRUE),
    mean_elapsed_total = mean(okdf$elapsed_total, na.rm = TRUE),
    median_elapsed_total = median(okdf$elapsed_total, na.rm = TRUE),
    mean_fval = mean(okdf$fval, na.rm = TRUE),
    median_fval = median(okdf$fval, na.rm = TRUE),
    mean_sig_main = mean(okdf$sig_main, na.rm = TRUE),
    mean_sig_aux = mean(okdf$sig_aux, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  cfg <- parse_args(args)
  suppressPackageStartupMessages(library(npRmpi))
  suppressPackageStartupMessages(library(microbenchmark))

  npRmpi.init(nslaves = cfg$nslaves)
  on.exit(try(npRmpi.quit(force = TRUE), silent = TRUE), add = TRUE)

  options(npRmpi.autodispatch = TRUE, np.messages = FALSE)

  seeds <- make_seeds(cfg)
  if (cfg$show_progress) {
    cat("Running one-off MPI benchmark\n")
    cat("fun=", cfg$fun,
        " n=", cfg$n,
        " nslaves=", cfg$nslaves,
        " times=", cfg$times,
        " seed_policy=", cfg$seed_policy,
        " real_data=", is_real_data_fun(cfg$fun), "\n", sep = "")
  }

  rows <- vector("list", length(seeds))
  for (i in seq_along(seeds)) {
    if (cfg$show_progress) cat("  run", i, "seed", seeds[[i]], "\n")
    rows[[i]] <- run_seed(cfg, seeds[[i]])
  }

  raw <- do.call(rbind, rows)
  summary <- summarize_results(raw)

  write.table(raw, file = cfg$out_raw, sep = ",", row.names = FALSE,
              col.names = !file.exists(cfg$out_raw), append = file.exists(cfg$out_raw),
              qmethod = "double")

  if (nrow(summary) > 0L) {
    write.table(summary, file = cfg$out_summary, sep = ",", row.names = FALSE,
                col.names = !file.exists(cfg$out_summary), append = file.exists(cfg$out_summary),
                qmethod = "double")
  }

  if (cfg$show_progress) {
    cat("Raw:", cfg$out_raw, "\n")
    cat("Summary:", cfg$out_summary, "\n")
    if (nrow(summary) > 0L) print(summary)
  }
}

if (sys.nframe() == 0) main()
