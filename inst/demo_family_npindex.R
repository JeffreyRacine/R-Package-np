npindex_demo_source_utils <- function() {
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

npindex_demo_env <- function(name, default = "") {
  value <- Sys.getenv(name, default)
  if (length(value) != 1L) value <- value[[1L]]
  value
}

npindex_demo_bool <- function(value, name) {
  value <- trimws(value)
  if (value %in% c("TRUE", "true", "1")) return(TRUE)
  if (value %in% c("FALSE", "false", "0")) return(FALSE)
  stop(name, " must be TRUE/FALSE", call. = FALSE)
}

npindex_demo_int <- function(value, name, required = TRUE) {
  value <- trimws(value)
  if (!nzchar(value)) {
    if (required) stop(name, " is required", call. = FALSE)
    return(NA_integer_)
  }
  out <- suppressWarnings(as.integer(value))
  if (is.na(out) || as.character(out) != value || out < 0L) {
    stop(name, " must be a non-negative integer", call. = FALSE)
  }
  out
}

npindex_demo_numeric_field <- function(x) {
  if (is.null(x) || !length(x)) return(NA_real_)
  suppressWarnings(as.numeric(x[[1L]]))
}

npindex_demo_matrix_path <- function() {
  explicit <- npindex_demo_env("NP_DEMO_MATRIX")
  if (nzchar(explicit)) return(explicit)
  tier <- npindex_demo_env("NP_DEMO_TIER", "smoke")
  src <- npindex_demo_env("NP_DEMO_SRC")
  candidates <- c(
    if (nzchar(src)) file.path(src, "..", "inst", "demo_matrices",
                               paste0("npindex-", tier, ".csv")) else "",
    system.file("demo_matrices", paste0("npindex-", tier, ".csv"),
                package = "npRmpi")
  )
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (!length(candidates)) {
    stop("Could not locate npindex demo matrix for tier=", tier, call. = FALSE)
  }
  candidates[[1L]]
}

npindex_demo_matrix <- function() {
  path <- npindex_demo_matrix_path()
  mat <- read.csv(path, stringsAsFactors = FALSE, na.strings = character())
  mat[is.na(mat)] <- ""
  required <- c("family", "case", "tier", "enabled", "default_n", "floor_n",
                "method", "regtype", "nomad", "degree", "degree.max",
                "bwtype", "gradients", "why")
  missing <- setdiff(required, names(mat))
  if (length(missing)) {
    stop("npindex demo matrix missing columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  mat <- mat[mat$enabled == "TRUE", , drop = FALSE]
  cases <- npindex_demo_env("NP_DEMO_CASES")
  if (!nzchar(cases)) cases <- npindex_demo_env("NP_DEMO_CASE")
  if (nzchar(cases)) {
    keep <- strsplit(cases, "[[:space:]]+", perl = TRUE)[[1L]]
    keep <- keep[nzchar(keep)]
    mat <- mat[mat$case %in% keep, , drop = FALSE]
  }
  if (!nrow(mat)) stop("npindex demo matrix selected no enabled rows", call. = FALSE)
  mat
}

npindex_demo_validate_row <- function(row) {
  if (!identical(row$family, "npindex")) {
    stop("npindex family script received family=", row$family, call. = FALSE)
  }
  if (!(row$method %in% c("ichimura", "kleinspady"))) {
    stop("unsupported npindex method=", row$method, call. = FALSE)
  }
  if (!(row$regtype %in% c("lc", "ll", "lp"))) {
    stop("unsupported npindex regtype=", row$regtype, call. = FALSE)
  }
  if (!identical(row$bwtype, "fixed")) {
    stop("npindex demo matrix currently supports bwtype=fixed only", call. = FALSE)
  }
  nomad <- npindex_demo_bool(row$nomad, "nomad")
  if (nomad && !identical(row$regtype, "lp")) {
    stop("npindex nomad rows must use regtype=lp", call. = FALSE)
  }
  if (nomad && nzchar(row$degree)) {
    stop("npindex nomad rows must leave degree empty for degree search",
         call. = FALSE)
  }
  if (nomad && !nzchar(row$degree.max)) {
    stop("npindex nomad rows must provide degree.max", call. = FALSE)
  }
  if (!nomad && identical(row$regtype, "lp") && !nzchar(row$degree)) {
    stop("npindex lp rows must provide degree when nomad=FALSE", call. = FALSE)
  }
  invisible(TRUE)
}

npindex_demo_data <- function(method, n) {
  set.seed(42)
  if (identical(method, "kleinspady")) {
    x <- rchisq(n, df = 3)
    x1 <- (ifelse(x < 6, x, 6) - 2.348) / 1.511
    x <- rnorm(n)
    x2 <- ifelse(abs(x) < 2, x, 2) / 0.8796
    y <- ifelse(x1 + x2 + rnorm(n) > 0, 1, 0)
  } else {
    x1 <- runif(n, min = -1, max = 1)
    x2 <- runif(n, min = -1, max = 1)
    y <- x1 - x2 + rnorm(n)
  }
  data.frame(x1 = x1, x2 = x2, y = y)
}

npindex_demo_run_row <- function(row, mode) {
  npindex_demo_validate_row(row)
  default_n <- npindex_demo_int(row$default_n, "default_n")
  floor_n <- npindex_demo_int(row$floor_n, "floor_n")
  nomad <- npindex_demo_bool(row$nomad, "nomad")
  gradients <- npindex_demo_bool(row$gradients, "gradients")
  degree <- if (identical(row$regtype, "lp") && !nomad) {
    npindex_demo_int(row$degree, "degree")
  } else {
    NA_integer_
  }
  degree.max <- if (nomad) npindex_demo_int(row$degree.max, "degree.max") else NA_integer_
  n <- np_demo_n(default_n, floor = floor_n)
  mydat <- npindex_demo_data(row$method, n)

  bw_args <- list(formula = y ~ x1 + x2, data = mydat, method = row$method,
                  regtype = row$regtype, bwtype = row$bwtype)
  if (identical(row$regtype, "lp") && !nomad) bw_args$degree <- degree
  if (nomad) {
    bw_args$nomad <- TRUE
    bw_args$degree.max <- degree.max
    bw_args$nmulti <- 1L
  }

  if (identical(mode, "profile")) {
    mpi.bcast.Robj2slave(mydat)
    mpi.bcast.Robj2slave(bw_args)
    mpi.bcast.Robj2slave(gradients)
    t <- system.time(mpi.bcast.cmd(bw <- do.call(npindexbw, bw_args),
                                    caller.execute = TRUE))
    t <- t + system.time(mpi.bcast.cmd(model <- npindex(bws = bw, gradients = gradients),
                                       caller.execute = TRUE))
  } else {
    t <- system.time(bw <- do.call(npindexbw, bw_args))
    t <- t + system.time(model <- npindex(bws = bw, gradients = gradients))
  }

  summary(bw)
  summary(model)
  cat("Elapsed time =", t[["elapsed"]], "\n")
  np_demo_result(row$case, mode, n, default_n, t[["elapsed"]],
                 family = row$family,
                 case = row$case,
                 tier = row$tier,
                 method = row$method,
                 regtype = row$regtype,
                 bwmethod = "NA",
                 nomad = nomad,
                 degree = if (is.na(degree)) NA else degree,
                 degree.max = if (is.na(degree.max)) NA else degree.max,
                 selected.degree = npindex_demo_numeric_field(bw$degree),
                 bwtype = row$bwtype,
                 gradients = gradients,
                 why = row$why,
                 fitted.length = length(model$mean),
                 objective = npindex_demo_numeric_field(bw$fval),
                 bandwidth = npindex_demo_numeric_field(bw$bw),
                 beta.second = if (length(bw$beta) >= 2L) as.numeric(bw$beta[[2L]]) else NA_real_,
                 num.feval = npindex_demo_numeric_field(bw$num.feval),
                 num.feval.fast = npindex_demo_numeric_field(bw$num.feval.fast),
                 nomad.time = npindex_demo_numeric_field(bw$nomad.time),
                 powell.time = npindex_demo_numeric_field(bw$powell.time))
}

npindex_demo_run_matrix <- function(mode) {
  npindex_demo_source_utils()
  mat <- npindex_demo_matrix()
  for (i in seq_len(nrow(mat))) {
    row <- mat[i, , drop = FALSE]
    cat("[npindex][", mode, "] ", row$case, "\n", sep = "")
    npindex_demo_run_row(row, mode)
  }
  invisible(TRUE)
}
