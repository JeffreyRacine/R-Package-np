npreg_demo_source_utils <- function() {
  if (exists("np_demo_n", mode = "function") &&
      exists("np_demo_result", mode = "function")) {
    return(invisible(TRUE))
  }
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

npreg_demo_env <- function(name, default = "") {
  value <- Sys.getenv(name, default)
  if (length(value) != 1L) value <- value[[1L]]
  value
}

npreg_demo_bool <- function(value, name) {
  value <- trimws(value)
  if (value %in% c("TRUE", "true", "1")) return(TRUE)
  if (value %in% c("FALSE", "false", "0")) return(FALSE)
  stop(name, " must be TRUE/FALSE", call. = FALSE)
}

npreg_demo_int <- function(value, name, required = TRUE) {
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

npreg_demo_numeric_field <- function(x) {
  if (is.null(x) || !length(x)) return(NA_real_)
  suppressWarnings(as.numeric(x[[1L]]))
}

npreg_demo_matrix_path <- function() {
  explicit <- npreg_demo_env("NP_DEMO_MATRIX")
  if (nzchar(explicit)) return(explicit)
  tier <- npreg_demo_env("NP_DEMO_TIER", "smoke")
  src <- npreg_demo_env("NP_DEMO_SRC")
  candidates <- c(
    if (nzchar(src)) file.path(src, "..", "inst", "demo_matrices",
                               paste0("npreg-", tier, ".csv")) else "",
    system.file("demo_matrices", paste0("npreg-", tier, ".csv"),
                package = "npRmpi")
  )
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (!length(candidates)) {
    stop("Could not locate npreg demo matrix for tier=", tier, call. = FALSE)
  }
  candidates[[1L]]
}

npreg_demo_matrix <- function() {
  path <- npreg_demo_matrix_path()
  mat <- read.csv(path, stringsAsFactors = FALSE, na.strings = character())
  mat[is.na(mat)] <- ""
  required <- c("family", "case", "tier", "enabled", "default_n", "floor_n",
                "regtype", "bwmethod", "nomad", "degree", "bwtype", "why")
  missing <- setdiff(required, names(mat))
  if (length(missing)) {
    stop("npreg demo matrix missing columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  mat <- mat[identical(TRUE, TRUE) & mat$enabled == "TRUE", , drop = FALSE]
  cases <- npreg_demo_env("NP_DEMO_CASES")
  if (!nzchar(cases)) cases <- npreg_demo_env("NP_DEMO_CASE")
  if (nzchar(cases)) {
    keep <- strsplit(cases, "[[:space:]]+", perl = TRUE)[[1L]]
    keep <- keep[nzchar(keep)]
    mat <- mat[mat$case %in% keep, , drop = FALSE]
  }
  if (!nrow(mat)) stop("npreg demo matrix selected no enabled rows", call. = FALSE)
  mat
}

npreg_demo_validate_row <- function(row) {
  if (!identical(row$family, "npreg")) {
    stop("npreg family script received family=", row$family, call. = FALSE)
  }
  if (!(row$regtype %in% c("lc", "ll", "lp"))) {
    stop("unsupported npreg regtype=", row$regtype, call. = FALSE)
  }
  if (!(row$bwmethod %in% c("cv.ls", "cv.aic"))) {
    stop("unsupported npreg bwmethod=", row$bwmethod, call. = FALSE)
  }
  if (!identical(row$bwtype, "fixed")) {
    stop("npreg demo matrix currently supports bwtype=fixed only", call. = FALSE)
  }
  nomad <- npreg_demo_bool(row$nomad, "nomad")
  if (nomad && !identical(row$regtype, "lp")) {
    stop("npreg nomad rows must use regtype=lp", call. = FALSE)
  }
  if (nomad && nzchar(row$degree)) {
    stop("npreg nomad rows must leave degree empty for degree search",
         call. = FALSE)
  }
  if (!nomad && identical(row$regtype, "lp") && !nzchar(row$degree)) {
    stop("npreg lp rows must provide degree when nomad=FALSE", call. = FALSE)
  }
  if (!nomad && !identical(row$regtype, "lp") && nzchar(row$degree)) {
    stop("npreg degree must be empty unless regtype=lp and nomad=FALSE",
         call. = FALSE)
  }
  invisible(TRUE)
}

npreg_demo_run_row <- function(row, mode) {
  npreg_demo_validate_row(row)
  default_n <- npreg_demo_int(row$default_n, "default_n")
  floor_n <- npreg_demo_int(row$floor_n, "floor_n")
  nomad <- npreg_demo_bool(row$nomad, "nomad")
  degree <- if (identical(row$regtype, "lp") && !nomad) {
    npreg_demo_int(row$degree, "degree")
  } else {
    NA_integer_
  }
  n <- np_demo_n(default_n, floor = floor_n)

  set.seed(42)
  x <- runif(n)
  z <- factor(rbinom(n, 1L, 0.5))
  y <- cos(2 * pi * x) + as.integer(z) - 1L + rnorm(n, sd = 0.25)
  mydat <- data.frame(y = y, x = x, z = z)
  rm(x, y, z)

  bw_args <- list(formula = y ~ x + z, data = mydat, regtype = row$regtype,
                  bwmethod = row$bwmethod, bwtype = row$bwtype)
  if (identical(row$regtype, "lp") && !nomad) bw_args$degree <- degree
  if (nomad) {
    bw_args$nomad <- TRUE
    bw_args$nmulti <- 1L
  }

  if (identical(mode, "profile")) {
    mpi.bcast.Robj2slave(mydat)
    mpi.bcast.Robj2slave(bw_args)
    t <- system.time(mpi.bcast.cmd(bw <- do.call(npregbw, bw_args),
                                    caller.execute = TRUE))
    t <- t + system.time(mpi.bcast.cmd(model <- npreg(bws = bw, data = mydat),
                                       caller.execute = TRUE))
  } else {
    t <- system.time(bw <- do.call(npregbw, bw_args))
    t <- t + system.time(model <- npreg(bws = bw, data = mydat))
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
                 selected.degree = npreg_demo_numeric_field(bw$degree),
                 bwtype = row$bwtype,
                 why = row$why,
                 fitted.length = length(model$mean),
                 objective = npreg_demo_numeric_field(bw$fval),
                 num.feval = npreg_demo_numeric_field(bw$num.feval),
                 num.feval.fast = npreg_demo_numeric_field(bw$num.feval.fast),
                 nomad.time = npreg_demo_numeric_field(bw$nomad.time),
                 powell.time = npreg_demo_numeric_field(bw$powell.time))
}

npreg_demo_run_matrix <- function(mode) {
  npreg_demo_source_utils()
  mat <- npreg_demo_matrix()
  for (i in seq_len(nrow(mat))) {
    row <- mat[i, , drop = FALSE]
    cat("[npreg][", mode, "] ", row$case, "\n", sep = "")
    npreg_demo_run_row(row, mode)
  }
  invisible(TRUE)
}
