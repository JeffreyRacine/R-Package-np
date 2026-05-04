npscoef_demo_source_utils <- function() {
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

npscoef_demo_env <- function(name, default = "") {
  value <- Sys.getenv(name, default)
  if (length(value) != 1L) value <- value[[1L]]
  value
}

npscoef_demo_bool <- function(value, name) {
  value <- trimws(value)
  if (value %in% c("TRUE", "true", "1")) return(TRUE)
  if (value %in% c("FALSE", "false", "0")) return(FALSE)
  stop(name, " must be TRUE/FALSE", call. = FALSE)
}

npscoef_demo_int <- function(value, name, required = TRUE) {
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

npscoef_demo_numeric_field <- function(...) {
  values <- list(...)
  for (value in values) {
    if (is.null(value) || !length(value)) next
    out <- suppressWarnings(as.numeric(value[[1L]]))
    if (!is.na(out)) return(out)
  }
  NA_real_
}

npscoef_demo_matrix_path <- function() {
  explicit <- npscoef_demo_env("NP_DEMO_MATRIX")
  if (nzchar(explicit)) return(explicit)
  tier <- npscoef_demo_env("NP_DEMO_TIER", "smoke")
  src <- npscoef_demo_env("NP_DEMO_SRC")
  candidates <- c(
    if (nzchar(src)) file.path(src, "..", "inst", "demo_matrices",
                               paste0("npscoef-", tier, ".csv")) else "",
    system.file("demo_matrices", paste0("npscoef-", tier, ".csv"),
                package = "npRmpi")
  )
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (!length(candidates)) {
    stop("Could not locate npscoef demo matrix for tier=", tier, call. = FALSE)
  }
  candidates[[1L]]
}

npscoef_demo_matrix <- function() {
  path <- npscoef_demo_matrix_path()
  mat <- read.csv(path, stringsAsFactors = FALSE, na.strings = character())
  mat[is.na(mat)] <- ""
  required <- c("family", "case", "tier", "enabled", "default_n", "floor_n",
                "regtype", "bwmethod", "nomad", "degree", "degree.max",
                "bwtype", "gradients", "why")
  missing <- setdiff(required, names(mat))
  if (length(missing)) {
    stop("npscoef demo matrix missing columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  mat <- mat[mat$enabled == "TRUE", , drop = FALSE]
  cases <- npscoef_demo_env("NP_DEMO_CASES")
  if (!nzchar(cases)) cases <- npscoef_demo_env("NP_DEMO_CASE")
  if (nzchar(cases)) {
    keep <- strsplit(cases, "[[:space:]]+", perl = TRUE)[[1L]]
    keep <- keep[nzchar(keep)]
    mat <- mat[mat$case %in% keep, , drop = FALSE]
  }
  if (!nrow(mat)) stop("npscoef demo matrix selected no enabled rows", call. = FALSE)
  mat
}

npscoef_demo_validate_row <- function(row) {
  if (!identical(row$family, "npscoef")) {
    stop("npscoef family script received family=", row$family, call. = FALSE)
  }
  if (!(row$regtype %in% c("lc", "ll", "lp"))) {
    stop("unsupported npscoef regtype=", row$regtype, call. = FALSE)
  }
  if (!identical(row$bwmethod, "cv.ls")) {
    stop("npscoef demo matrix currently supports bwmethod=cv.ls only",
         call. = FALSE)
  }
  if (!identical(row$bwtype, "fixed")) {
    stop("npscoef demo matrix currently supports bwtype=fixed only",
         call. = FALSE)
  }
  nomad <- npscoef_demo_bool(row$nomad, "nomad")
  if (nomad && !identical(row$regtype, "lp")) {
    stop("npscoef nomad rows must use regtype=lp", call. = FALSE)
  }
  if (nomad && nzchar(row$degree)) {
    stop("npscoef nomad rows must leave degree empty for degree search",
         call. = FALSE)
  }
  if (nomad && !nzchar(row$degree.max)) {
    stop("npscoef nomad rows must provide degree.max", call. = FALSE)
  }
  if (!nomad && identical(row$regtype, "lp") && !nzchar(row$degree)) {
    stop("npscoef lp rows must provide degree when nomad=FALSE", call. = FALSE)
  }
  invisible(TRUE)
}

npscoef_demo_data <- function(n) {
  set.seed(42)
  x <- runif(n)
  z <- runif(n, min = -2, max = 2)
  y <- x * exp(z) * (1.0 + rnorm(n, sd = 0.2))
  data.frame(x = x, y = y, z = z)
}

npscoef_demo_run_row <- function(row, mode) {
  npscoef_demo_validate_row(row)
  default_n <- npscoef_demo_int(row$default_n, "default_n")
  floor_n <- npscoef_demo_int(row$floor_n, "floor_n")
  nomad <- npscoef_demo_bool(row$nomad, "nomad")
  gradients <- npscoef_demo_bool(row$gradients, "gradients")
  degree <- if (identical(row$regtype, "lp") && !nomad) {
    npscoef_demo_int(row$degree, "degree")
  } else {
    NA_integer_
  }
  degree.max <- if (nomad) npscoef_demo_int(row$degree.max, "degree.max") else NA_integer_
  n <- np_demo_n(default_n, floor = floor_n)
  mydat <- npscoef_demo_data(n)

  bw_args <- list(formula = y ~ x | z, data = mydat, regtype = row$regtype,
                  bwmethod = row$bwmethod, bwtype = row$bwtype)
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
    t <- system.time(mpi.bcast.cmd(bw <- do.call(npscoefbw, bw_args),
                                   caller.execute = TRUE))
    t <- t + system.time(mpi.bcast.cmd(model <- npscoef(bws = bw, gradients = gradients),
                                       caller.execute = TRUE))
  } else {
    t <- system.time(bw <- do.call(npscoefbw, bw_args))
    t <- t + system.time(model <- npscoef(bws = bw, gradients = gradients))
  }

  summary(bw)
  summary(model)
  cat("Elapsed time =", t[["elapsed"]], "\n")
  np_demo_result(row$case, mode, n, default_n, t[["elapsed"]],
                 family = row$family,
                 case = row$case,
                 tier = row$tier,
                 regtype = row$regtype,
                 bwmethod = row$bwmethod,
                 nomad = nomad,
                 degree = if (is.na(degree)) NA else degree,
                 degree.max = if (is.na(degree.max)) NA else degree.max,
                 selected.degree = npscoef_demo_numeric_field(bw$degree),
                 bwtype = row$bwtype,
                 gradients = gradients,
                 why = row$why,
                 fitted.length = length(model$mean),
                 objective = npscoef_demo_numeric_field(bw$fval),
                 bandwidth = npscoef_demo_numeric_field(bw$bw),
                 num.feval = npscoef_demo_numeric_field(bw$num.feval,
                                                        bw$num.feval.function),
                 num.feval.fast = npscoef_demo_numeric_field(bw$num.feval.fast),
                 nomad.time = npscoef_demo_numeric_field(bw$nomad.time),
                 powell.time = npscoef_demo_numeric_field(bw$powell.time))
}

npscoef_demo_run_matrix <- function(mode) {
  npscoef_demo_source_utils()
  mat <- npscoef_demo_matrix()
  for (i in seq_len(nrow(mat))) {
    row <- mat[i, , drop = FALSE]
    cat("[npscoef][", mode, "] ", row$case, "\n", sep = "")
    npscoef_demo_run_row(row, mode)
  }
  invisible(TRUE)
}
