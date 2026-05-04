npaux_demo_source_utils <- function() {
  if (exists("np_demo_n", mode = "function") &&
      exists("np_demo_result", mode = "function")) return(invisible(TRUE))
  src <- Sys.getenv("NP_DEMO_SRC", "")
  candidates <- c(
    Sys.getenv("NP_DEMO_UTILS", ""),
    if (nzchar(src)) file.path(src, "..", "inst", "demo_utils.R") else "",
    if (nzchar(src)) file.path(src, "demo_utils.R") else "",
    system.file("demo_utils.R", package = "npRmpi")
  )
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (!length(candidates)) stop("Could not locate demo_utils.R", call. = FALSE)
  source(candidates[[1L]])
  invisible(TRUE)
}

npaux_demo_env <- function(name, default = "") {
  value <- Sys.getenv(name, default)
  if (length(value) != 1L) value <- value[[1L]]
  value
}

npaux_demo_int <- function(value, name) {
  value <- trimws(value)
  out <- suppressWarnings(as.integer(value))
  if (!nzchar(value) || is.na(out) || as.character(out) != value || out <= 0L) {
    stop(name, " must be a positive integer", call. = FALSE)
  }
  out
}

npaux_demo_num <- function(x, i = 1L) {
  raw <- tryCatch(unlist(x, recursive = TRUE, use.names = FALSE),
                  error = function(e) numeric(0))
  out <- tryCatch(suppressWarnings(as.numeric(raw)),
                  error = function(e) numeric(0))
  out <- out[is.finite(out)]
  if (length(out) < i) NA_real_ else out[[i]]
}

npaux_demo_matrix_path <- function(family) {
  explicit <- npaux_demo_env("NP_DEMO_MATRIX")
  if (nzchar(explicit)) return(explicit)
  tier <- npaux_demo_env("NP_DEMO_TIER", "smoke")
  src <- npaux_demo_env("NP_DEMO_SRC")
  candidates <- c(
    if (nzchar(src)) file.path(src, "..", "inst", "demo_matrices",
                               paste0(family, "-", tier, ".csv")) else "",
    system.file("demo_matrices", paste0(family, "-", tier, ".csv"),
                package = "npRmpi")
  )
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (!length(candidates)) {
    stop("Could not locate ", family, " demo matrix for tier=", tier,
         call. = FALSE)
  }
  candidates[[1L]]
}

npaux_demo_matrix <- function(family) {
  path <- npaux_demo_matrix_path(family)
  mat <- read.csv(path, stringsAsFactors = FALSE, na.strings = character())
  mat[is.na(mat)] <- ""
  required <- c("family", "case", "tier", "enabled", "default_n", "floor_n",
                "why")
  missing <- setdiff(required, names(mat))
  if (length(missing)) {
    stop(family, " demo matrix missing columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  mat <- mat[mat$enabled == "TRUE", , drop = FALSE]
  cases <- npaux_demo_env("NP_DEMO_CASES")
  if (!nzchar(cases)) cases <- npaux_demo_env("NP_DEMO_CASE")
  if (nzchar(cases)) {
    keep <- strsplit(cases, "[[:space:]]+", perl = TRUE)[[1L]]
    keep <- keep[nzchar(keep)]
    mat <- mat[mat$case %in% keep, , drop = FALSE]
  }
  if (!nrow(mat)) stop(family, " demo matrix selected no enabled rows", call. = FALSE)
  mat
}

npaux_demo_payload <- function(family, row, n) {
  set.seed(42)
  switch(
    family,
    npcmstest = {
      x1 <- rnorm(n)
      x2 <- runif(n)
      y <- 1.0 + x1 + x2 + rnorm(n)
      X <- data.frame(x1 = x1, x2 = x2)
      list(model = lm(y ~ x1 + x2, x = TRUE, y = TRUE),
           X = X, y = y,
           boot.num = npaux_demo_int(row$boot.num, "boot.num"))
    },
    npconmode = {
      x1 <- factor(rbinom(n, 1, .5))
      x2 <- ordered(sample(0:2, n, replace = TRUE))
      x3 <- rnorm(n)
      low <- factor(ifelse(x3 + as.numeric(x1) + rnorm(n) > 0, 1, 0))
      list(data = data.frame(low = low, x1 = x1, x2 = x2, x3 = x3))
    },
    npcopula = {
      rho <- 0.8
      x <- rnorm(n)
      y <- rho * x + sqrt(1.0 - rho^2) * rnorm(n)
      list(data = data.frame(x = x, y = y),
           grid = data.frame(x = seq(0, 1, length.out = 10L),
                             y = seq(0, 1, length.out = 10L)))
    },
    npqreg = {
      x <- rnorm(n)
      y <- x + rnorm(n)
      list(data = data.frame(y = y, x = x))
    },
    npregiv = {
      v <- rnorm(n, mean = 0, sd = 0.27)
      eps <- rnorm(n, mean = 0, sd = 0.05)
      u <- -0.5 * v + eps
      w <- rnorm(n)
      z <- 0.2 * w + v
      y <- z^2 + u
      list(y = y, z = z, w = w)
    },
    stop("unsupported auxiliary demo family=", family, call. = FALSE)
  )
}

npaux_demo_execute <- function(family, payload) {
  switch(
    family,
    npcmstest = {
      output <- npcmstest(model = payload$model, xdat = payload$X,
                          ydat = payload$y, nmulti = 1, ftol = .01, tol = .01,
                          boot.num = payload$boot.num)
      list(output = output,
           estimate.first = npaux_demo_num(output, 1L),
           estimate.second = npaux_demo_num(output, 2L),
           fitted.length = length(payload$y))
    },
    npconmode = {
      bw <- npcdensbw(low ~ x1 + x2 + x3, data = payload$data,
                      bwmethod = "cv.ls")
      fit <- npconmode(bws = bw)
      list(output = fit, bw = bw,
           objective = npaux_demo_num(bw$fval, 1L),
           bandwidth = {
             out <- npaux_demo_num(bw$bw, 1L)
             if (is.na(out)) npaux_demo_num(bw$bandwidth, 1L) else out
           },
           estimate.first = npaux_demo_num(fit$conmode, 1L),
           estimate.second = npaux_demo_num(fit$condens, 1L),
           fitted.length = length(fit$conmode))
    },
    npcopula = {
      bw <- npudistbw(~x + y, data = payload$data)
      fit <- npcopula(bws = bw, data = payload$data, u = payload$grid)
      list(output = fit, bw = bw,
           objective = npaux_demo_num(bw$fval, 1L),
           bandwidth = {
             out <- npaux_demo_num(bw$bw, 1L)
             if (is.na(out)) npaux_demo_num(bw$bandwidth, 1L) else out
           },
           estimate.first = npaux_demo_num(fit, 1L),
           estimate.second = npaux_demo_num(fit, 2L),
           fitted.length = NROW(fit))
    },
    npqreg = {
      bw <- npcdistbw(y ~ x, data = payload$data)
      q25 <- npqreg(bws = bw, tau = 0.25)
      q50 <- npqreg(bws = bw, tau = 0.50)
      q75 <- npqreg(bws = bw, tau = 0.75)
      list(output = q50, bw = bw,
           objective = npaux_demo_num(bw$fval, 1L),
           bandwidth = {
             out <- npaux_demo_num(bw$bw, 1L)
             if (is.na(out)) npaux_demo_num(bw$bandwidth, 1L) else out
           },
           estimate.first = npaux_demo_num(q50$quantile, 1L),
           estimate.second = npaux_demo_num(q25$quantile, 1L),
           fitted.length = length(q75$quantile))
    },
    npregiv = {
      fit <- npregiv(y = payload$y, z = payload$z, w = payload$w)
      list(output = fit,
           estimate.first = npaux_demo_num(fit$phi, 1L),
           estimate.second = npaux_demo_num(fit$phi, 2L),
           fitted.length = length(fit$phi))
    },
    stop("unsupported auxiliary demo family=", family, call. = FALSE)
  )
}

npaux_demo_run_row <- function(family, row, mode) {
  if (!identical(row$family, family)) {
    stop("auxiliary family script received family=", row$family,
         " while running ", family, call. = FALSE)
  }
  default_n <- npaux_demo_int(row$default_n, "default_n")
  floor_n <- npaux_demo_int(row$floor_n, "floor_n")
  n <- np_demo_n(default_n, floor = floor_n)
  payload <- npaux_demo_payload(family, row, n)

  if (identical(mode, "profile")) {
    mpi.bcast.Robj2slave(family)
    mpi.bcast.Robj2slave(payload)
    t <- system.time(mpi.bcast.cmd(result <- get("npaux_demo_execute", envir = .GlobalEnv)(family, payload),
                                   caller.execute = TRUE))
  } else {
    t <- system.time(result <- npaux_demo_execute(family, payload))
  }

  if (!is.null(result$bw)) summary(result$bw)
  print(result$output)
  cat("Elapsed time =", t[["elapsed"]], "\n")
  np_demo_result(row$case, mode, n, default_n, t[["elapsed"]],
                 family = family,
                 case = row$case,
                 tier = row$tier,
                 why = row$why,
                 fitted.length = result$fitted.length,
                 objective = if (is.null(result$objective)) NA_real_ else result$objective,
                 bandwidth = if (is.null(result$bandwidth)) NA_real_ else result$bandwidth,
                 estimate.first = result$estimate.first,
                 estimate.second = result$estimate.second)
}

npaux_demo_run_matrix <- function(family, mode) {
  npaux_demo_source_utils()
  mat <- npaux_demo_matrix(family)
  for (i in seq_len(nrow(mat))) {
    row <- mat[i, , drop = FALSE]
    cat("[", family, "][", mode, "] ", row$case, "\n", sep = "")
    npaux_demo_run_row(family, row, mode)
  }
  invisible(TRUE)
}
