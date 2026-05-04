npplreg_demo_source_utils <- function() {
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

npplreg_demo_env <- function(name, default = "") {
  value <- Sys.getenv(name, default)
  if (length(value) != 1L) value <- value[[1L]]
  value
}

npplreg_demo_bool <- function(value, name) {
  value <- trimws(value)
  if (value %in% c("TRUE", "true", "1")) return(TRUE)
  if (value %in% c("FALSE", "false", "0")) return(FALSE)
  stop(name, " must be TRUE/FALSE", call. = FALSE)
}

npplreg_demo_int <- function(value, name, required = TRUE) {
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

npplreg_demo_numeric_field <- function(...) {
  values <- list(...)
  for (value in values) {
    if (is.null(value) || !length(value)) next
    out <- suppressWarnings(as.numeric(value[[1L]]))
    if (!is.na(out)) return(out)
  }
  NA_real_
}

npplreg_demo_numeric_vector <- function(value) {
  out <- suppressWarnings(as.numeric(unlist(value, use.names = FALSE)))
  out <- out[is.finite(out)]
  if (!length(out)) numeric(0) else out
}

npplreg_demo_matrix_path <- function() {
  explicit <- npplreg_demo_env("NP_DEMO_MATRIX")
  if (nzchar(explicit)) return(explicit)
  tier <- npplreg_demo_env("NP_DEMO_TIER", "smoke")
  src <- npplreg_demo_env("NP_DEMO_SRC")
  candidates <- c(
    if (nzchar(src)) file.path(src, "..", "inst", "demo_matrices",
                               paste0("npplreg-", tier, ".csv")) else "",
    system.file("demo_matrices", paste0("npplreg-", tier, ".csv"),
                package = "npRmpi")
  )
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (!length(candidates)) {
    stop("Could not locate npplreg demo matrix for tier=", tier, call. = FALSE)
  }
  candidates[[1L]]
}

npplreg_demo_matrix <- function() {
  path <- npplreg_demo_matrix_path()
  mat <- read.csv(path, stringsAsFactors = FALSE, na.strings = character())
  mat[is.na(mat)] <- ""
  required <- c("family", "case", "tier", "enabled", "default_n", "floor_n",
                "regtype", "bwmethod", "nomad", "degree", "degree.max",
                "bwtype", "why")
  missing <- setdiff(required, names(mat))
  if (length(missing)) {
    stop("npplreg demo matrix missing columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  mat <- mat[mat$enabled == "TRUE", , drop = FALSE]
  cases <- npplreg_demo_env("NP_DEMO_CASES")
  if (!nzchar(cases)) cases <- npplreg_demo_env("NP_DEMO_CASE")
  if (nzchar(cases)) {
    keep <- strsplit(cases, "[[:space:]]+", perl = TRUE)[[1L]]
    keep <- keep[nzchar(keep)]
    mat <- mat[mat$case %in% keep, , drop = FALSE]
  }
  if (!nrow(mat)) stop("npplreg demo matrix selected no enabled rows", call. = FALSE)
  mat
}

npplreg_demo_validate_row <- function(row) {
  if (!identical(row$family, "npplreg")) {
    stop("npplreg family script received family=", row$family, call. = FALSE)
  }
  if (!(row$regtype %in% c("lc", "ll", "lp"))) {
    stop("unsupported npplreg regtype=", row$regtype, call. = FALSE)
  }
  if (!identical(row$bwmethod, "cv.ls")) {
    stop("npplreg demo matrix currently supports bwmethod=cv.ls only",
         call. = FALSE)
  }
  if (!identical(row$bwtype, "fixed")) {
    stop("npplreg demo matrix currently supports bwtype=fixed only",
         call. = FALSE)
  }
  nomad <- npplreg_demo_bool(row$nomad, "nomad")
  if (nomad && !identical(row$regtype, "lp")) {
    stop("npplreg nomad rows must use regtype=lp", call. = FALSE)
  }
  if (nomad && nzchar(row$degree)) {
    stop("npplreg nomad rows must leave degree empty for degree search",
         call. = FALSE)
  }
  if (nomad && !nzchar(row$degree.max)) {
    stop("npplreg nomad rows must provide degree.max", call. = FALSE)
  }
  if (!nomad && identical(row$regtype, "lp") && !nzchar(row$degree)) {
    stop("npplreg lp rows must provide degree when nomad=FALSE", call. = FALSE)
  }
  invisible(TRUE)
}

npplreg_demo_data <- function(n) {
  set.seed(42)
  x1 <- rnorm(n)
  x2 <- ordered(rbinom(n, 5, .3))
  z1 <- ordered(rbinom(n, 2, .3))
  z2 <- rnorm(n)
  y <- 1 + x1 + as.numeric(x2) + as.numeric(z1) + sin(z2) + rnorm(n)
  data.frame(y = y, x1 = x1, x2 = x2, z1 = z1, z2 = z2)
}

npplreg_demo_run_row <- function(row, mode) {
  npplreg_demo_validate_row(row)
  default_n <- npplreg_demo_int(row$default_n, "default_n")
  floor_n <- npplreg_demo_int(row$floor_n, "floor_n")
  nomad <- npplreg_demo_bool(row$nomad, "nomad")
  degree <- if (identical(row$regtype, "lp") && !nomad) {
    npplreg_demo_int(row$degree, "degree")
  } else {
    NA_integer_
  }
  degree.max <- if (nomad) npplreg_demo_int(row$degree.max, "degree.max") else NA_integer_
  n <- np_demo_n(default_n, floor = floor_n)
  mydat <- npplreg_demo_data(n)

  bw_args <- list(formula = y ~ x1 + x2 | z1 + z2, data = mydat,
                  regtype = row$regtype, bwmethod = row$bwmethod,
                  bwtype = row$bwtype)
  if (identical(row$regtype, "lp") && !nomad) bw_args$degree <- degree
  if (nomad) {
    bw_args$nomad <- TRUE
    bw_args$degree.max <- degree.max
    bw_args$nmulti <- 1L
  }

  if (identical(mode, "profile")) {
    mpi.bcast.Robj2slave(mydat)
    mpi.bcast.Robj2slave(bw_args)
    t <- system.time(mpi.bcast.cmd(bw <- do.call(npplregbw, bw_args),
                                   caller.execute = TRUE))
    t <- t + system.time(mpi.bcast.cmd(model <- npplreg(bws = bw),
                                       caller.execute = TRUE))
  } else {
    t <- system.time(bw <- do.call(npplregbw, bw_args))
    t <- t + system.time(model <- npplreg(bws = bw))
  }

  bandwidth <- npplreg_demo_numeric_vector(bw$bandwidth)
  coef <- npplreg_demo_numeric_vector(model$xcoef)
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
                 selected.degree = npplreg_demo_numeric_field(bw$degree),
                 bwtype = row$bwtype,
                 why = row$why,
                 fitted.length = length(model$mean),
                 objective = npplreg_demo_numeric_field(bw$fval),
                 bandwidth = if (length(bandwidth) >= 1L) bandwidth[[1L]] else NA_real_,
                 bandwidth.second = if (length(bandwidth) >= 2L) bandwidth[[2L]] else NA_real_,
                 coef.first = if (length(coef) >= 1L) coef[[1L]] else NA_real_,
                 coef.second = if (length(coef) >= 2L) coef[[2L]] else NA_real_,
                 num.feval = npplreg_demo_numeric_field(bw$num.feval),
                 num.feval.fast = npplreg_demo_numeric_field(bw$num.feval.fast),
                 nomad.time = npplreg_demo_numeric_field(bw$nomad.time),
                 powell.time = npplreg_demo_numeric_field(bw$powell.time))
}

npplreg_demo_run_matrix <- function(mode) {
  npplreg_demo_source_utils()
  mat <- npplreg_demo_matrix()
  for (i in seq_len(nrow(mat))) {
    row <- mat[i, , drop = FALSE]
    cat("[npplreg][", mode, "] ", row$case, "\n", sep = "")
    npplreg_demo_run_row(row, mode)
  }
  invisible(TRUE)
}
