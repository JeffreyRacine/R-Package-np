npudist_demo_source_utils <- function() {
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

npudist_demo_env <- function(name, default = "") {
  value <- Sys.getenv(name, default)
  if (length(value) != 1L) value <- value[[1L]]
  value
}

npudist_demo_bool <- function(value, name) {
  value <- trimws(value)
  if (value %in% c("TRUE", "true", "1")) return(TRUE)
  if (value %in% c("FALSE", "false", "0")) return(FALSE)
  stop(name, " must be TRUE/FALSE", call. = FALSE)
}

npudist_demo_int <- function(value, name) {
  value <- trimws(value)
  out <- suppressWarnings(as.integer(value))
  if (!nzchar(value) || is.na(out) || as.character(out) != value || out <= 0L) {
    stop(name, " must be a positive integer", call. = FALSE)
  }
  out
}

npudist_demo_numeric_field <- function(...) {
  values <- list(...)
  for (value in values) {
    if (is.null(value) || !length(value)) next
    out <- suppressWarnings(as.numeric(value[[1L]]))
    if (!is.na(out)) return(out)
  }
  NA_real_
}

npudist_demo_matrix_path <- function() {
  explicit <- npudist_demo_env("NP_DEMO_MATRIX")
  if (nzchar(explicit)) return(explicit)
  tier <- npudist_demo_env("NP_DEMO_TIER", "smoke")
  src <- npudist_demo_env("NP_DEMO_SRC")
  candidates <- c(
    if (nzchar(src)) file.path(src, "..", "inst", "demo_matrices",
                               paste0("npudist-", tier, ".csv")) else "",
    system.file("demo_matrices", paste0("npudist-", tier, ".csv"),
                package = "npRmpi")
  )
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (!length(candidates)) {
    stop("Could not locate npudist demo matrix for tier=", tier, call. = FALSE)
  }
  candidates[[1L]]
}

npudist_demo_matrix <- function() {
  path <- npudist_demo_matrix_path()
  mat <- read.csv(path, stringsAsFactors = FALSE, na.strings = character())
  mat[is.na(mat)] <- ""
  required <- c("family", "case", "tier", "enabled", "default_n", "floor_n",
                "bwmethod", "bwtype", "do.full.integral", "why")
  missing <- setdiff(required, names(mat))
  if (length(missing)) {
    stop("npudist demo matrix missing columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  mat <- mat[mat$enabled == "TRUE", , drop = FALSE]
  cases <- npudist_demo_env("NP_DEMO_CASES")
  if (!nzchar(cases)) cases <- npudist_demo_env("NP_DEMO_CASE")
  if (nzchar(cases)) {
    keep <- strsplit(cases, "[[:space:]]+", perl = TRUE)[[1L]]
    keep <- keep[nzchar(keep)]
    mat <- mat[mat$case %in% keep, , drop = FALSE]
  }
  if (!nrow(mat)) stop("npudist demo matrix selected no enabled rows", call. = FALSE)
  mat
}

npudist_demo_validate_row <- function(row) {
  if (!identical(row$family, "npudist")) {
    stop("npudist family script received family=", row$family, call. = FALSE)
  }
  if (!identical(row$bwmethod, "cv.cdf")) {
    stop("unsupported npudist bwmethod=", row$bwmethod, call. = FALSE)
  }
  if (!identical(row$bwtype, "fixed")) {
    stop("npudist demo matrix currently supports bwtype=fixed only",
         call. = FALSE)
  }
  invisible(TRUE)
}

npudist_demo_data <- function(n) {
  set.seed(42)
  data.frame(x = rnorm(n))
}

npudist_demo_run_row <- function(row, mode) {
  npudist_demo_validate_row(row)
  default_n <- npudist_demo_int(row$default_n, "default_n")
  floor_n <- npudist_demo_int(row$floor_n, "floor_n")
  do.full.integral <- npudist_demo_bool(row$do.full.integral, "do.full.integral")
  n <- np_demo_n(default_n, floor = floor_n)
  mydat <- npudist_demo_data(n)
  bw_args <- list(formula = ~x, data = mydat, bwmethod = row$bwmethod,
                  bwtype = row$bwtype, do.full.integral = do.full.integral)

  if (identical(mode, "profile")) {
    mpi.bcast.Robj2slave(mydat)
    mpi.bcast.Robj2slave(bw_args)
    t <- system.time(mpi.bcast.cmd(bw <- do.call(npudistbw, bw_args),
                                   caller.execute = TRUE))
    t <- t + system.time(mpi.bcast.cmd(fit <- npudist(bws = bw),
                                       caller.execute = TRUE))
  } else {
    t <- system.time(bw <- do.call(npudistbw, bw_args))
    t <- t + system.time(fit <- npudist(bws = bw))
  }

  summary(bw)
  summary(fit)
  cat("Elapsed time =", t[["elapsed"]], "\n")
  np_demo_result(row$case, mode, n, default_n, t[["elapsed"]],
                 family = row$family,
                 case = row$case,
                 tier = row$tier,
                 bwmethod = row$bwmethod,
                 bwtype = row$bwtype,
                 do.full.integral = do.full.integral,
                 why = row$why,
                 fitted.length = length(fit$dist),
                 objective = npudist_demo_numeric_field(bw$fval),
                 bandwidth = npudist_demo_numeric_field(bw$bw),
                 distribution.first = npudist_demo_numeric_field(fit$dist),
                 num.feval = npudist_demo_numeric_field(bw$num.feval),
                 num.feval.fast = npudist_demo_numeric_field(bw$num.feval.fast))
}

npudist_demo_run_matrix <- function(mode) {
  npudist_demo_source_utils()
  mat <- npudist_demo_matrix()
  for (i in seq_len(nrow(mat))) {
    row <- mat[i, , drop = FALSE]
    cat("[npudist][", mode, "] ", row$case, "\n", sep = "")
    npudist_demo_run_row(row, mode)
  }
  invisible(TRUE)
}
