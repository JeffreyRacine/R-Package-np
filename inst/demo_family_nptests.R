nptest_demo_source_utils <- function() {
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

nptest_demo_env <- function(name, default = "") {
  value <- Sys.getenv(name, default)
  if (length(value) != 1L) value <- value[[1L]]
  value
}

nptest_demo_bool <- function(value, name) {
  value <- trimws(value)
  if (value %in% c("TRUE", "true", "1")) return(TRUE)
  if (value %in% c("FALSE", "false", "0")) return(FALSE)
  stop(name, " must be TRUE/FALSE", call. = FALSE)
}

nptest_demo_int <- function(value, name) {
  value <- trimws(value)
  out <- suppressWarnings(as.integer(value))
  if (!nzchar(value) || is.na(out) || as.character(out) != value || out <= 0L) {
    stop(name, " must be a positive integer", call. = FALSE)
  }
  out
}

nptest_demo_numeric_vector <- function(x) {
  raw <- tryCatch(unlist(x, recursive = TRUE, use.names = FALSE),
                  error = function(e) numeric(0))
  out <- tryCatch(suppressWarnings(as.numeric(raw)),
                  error = function(e) numeric(0))
  out[is.finite(out)]
}

nptest_demo_numeric_field <- function(x, i = 1L) {
  x <- nptest_demo_numeric_vector(x)
  if (length(x) < i) NA_real_ else x[[i]]
}

nptest_demo_matrix_path <- function(family) {
  explicit <- nptest_demo_env("NP_DEMO_MATRIX")
  if (nzchar(explicit)) return(explicit)
  tier <- nptest_demo_env("NP_DEMO_TIER", "smoke")
  src <- nptest_demo_env("NP_DEMO_SRC")
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

nptest_demo_matrix <- function(family) {
  path <- nptest_demo_matrix_path(family)
  mat <- read.csv(path, stringsAsFactors = FALSE, na.strings = character())
  mat[is.na(mat)] <- ""
  required <- c("family", "case", "tier", "enabled", "default_n", "floor_n",
                "boot.num", "method", "why")
  missing <- setdiff(required, names(mat))
  if (length(missing)) {
    stop(family, " demo matrix missing columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  mat <- mat[mat$enabled == "TRUE", , drop = FALSE]
  cases <- nptest_demo_env("NP_DEMO_CASES")
  if (!nzchar(cases)) cases <- nptest_demo_env("NP_DEMO_CASE")
  if (nzchar(cases)) {
    keep <- strsplit(cases, "[[:space:]]+", perl = TRUE)[[1L]]
    keep <- keep[nzchar(keep)]
    mat <- mat[mat$case %in% keep, , drop = FALSE]
  }
  if (!nrow(mat)) stop(family, " demo matrix selected no enabled rows", call. = FALSE)
  mat
}

nptest_demo_ar_series <- function(phi, epsilon) {
  m <- length(epsilon)
  series <- numeric(m)
  series[1L] <- epsilon[1L] / (1.0 - phi)
  for (i in 2:m)
    series[i] <- phi * series[i - 1L] + epsilon[i]
  series
}

nptest_demo_payload <- function(family, row, n) {
  boot.num <- nptest_demo_int(row$boot.num, "boot.num")
  set.seed(42)
  switch(
    family,
    npdeneqtest = {
      list(fun = "npdeneqtest",
           args = list(x = data.frame(x = rnorm(n)),
                       y = data.frame(x = rnorm(n)),
                       boot.num = boot.num))
    },
    npdeptest = {
      x <- rnorm(n)
      y <- 1.0 + x + rnorm(n)
      list(fun = "npdeptest",
           args = list(data.x = y, data.y = fitted(lm(y ~ x)),
                       boot.num = boot.num, method = row$method))
    },
    npsdeptest = {
      list(fun = "npsdeptest",
           args = list(data = nptest_demo_ar_series(0.95, rnorm(n)),
                       lag.num = 2L, boot.num = boot.num, method = row$method))
    },
    npsymtest = {
      list(fun = "npsymtest",
           args = list(data = nptest_demo_ar_series(0.5, rnorm(n)),
                       boot.num = boot.num, boot.method = "geom",
                       method = row$method))
    },
    npunitest = {
      list(fun = "npunitest",
           args = list(data.x = rnorm(n), data.y = rnorm(n),
                       method = row$method, bootstrap = TRUE,
                       boot.num = boot.num))
    },
    npsigtest = {
      z <- factor(rbinom(n, 1, .5))
      x1 <- rnorm(n)
      x2 <- runif(n, -2, 2)
      y <- x1 + x2 + rnorm(n)
      list(fun = "npsigtest",
           args = list(mydat = data.frame(z = z, x1 = x1, x2 = x2, y = y),
                       boot.num = boot.num))
    },
    stop("unsupported test demo family=", family, call. = FALSE)
  )
}

nptest_demo_run_payload <- function(payload) {
  if (!identical(payload$fun, "npsigtest"))
    return(do.call(payload$fun, payload$args))
  mydat <- payload$args$mydat
  model <- npreg(y ~ z + x1 + x2, regtype = "ll", bwmethod = "cv.aic",
                 data = mydat)
  npsigtest(model, boot.num = payload$args$boot.num)
}

nptest_demo_run_row <- function(family, row, mode) {
  if (!identical(row$family, family)) {
    stop("test family script received family=", row$family,
         " while running ", family, call. = FALSE)
  }
  default_n <- nptest_demo_int(row$default_n, "default_n")
  floor_n <- nptest_demo_int(row$floor_n, "floor_n")
  boot.num <- nptest_demo_int(row$boot.num, "boot.num")
  n <- np_demo_n(default_n, floor = floor_n)
  payload <- nptest_demo_payload(family, row, n)

  if (identical(mode, "profile")) {
    mpi.bcast.Robj2slave(payload)
    t <- system.time(mpi.bcast.cmd(output <- if (identical(payload$fun, "npsigtest")) {
      mydat <- payload$args$mydat
      model <- npreg(y ~ z + x1 + x2, regtype = "ll", bwmethod = "cv.aic",
                     data = mydat)
      npsigtest(model, boot.num = payload$args$boot.num)
    } else {
      do.call(payload$fun, payload$args)
    }, caller.execute = TRUE))
  } else {
    t <- system.time(output <- nptest_demo_run_payload(payload))
  }

  print(output)
  cat("Elapsed time =", t[["elapsed"]], "\n")
  p.value <- if (!is.null(output$p.value)) output$p.value else
    if (!is.null(output$pval)) output$pval else NA_real_
  np_demo_result(row$case, mode, n, default_n, t[["elapsed"]],
                 family = family,
                 case = row$case,
                 tier = row$tier,
                 method = row$method,
                 boot.num = boot.num,
                 why = row$why,
                 statistic.first = nptest_demo_numeric_field(output, 1L),
                 statistic.second = nptest_demo_numeric_field(output, 2L),
                 p.value = nptest_demo_numeric_field(p.value, 1L))
}

nptest_demo_run_matrix <- function(family, mode) {
  nptest_demo_source_utils()
  mat <- nptest_demo_matrix(family)
  for (i in seq_len(nrow(mat))) {
    row <- mat[i, , drop = FALSE]
    cat("[", family, "][", mode, "] ", row$case, "\n", sep = "")
    nptest_demo_run_row(family, row, mode)
  }
  invisible(TRUE)
}
