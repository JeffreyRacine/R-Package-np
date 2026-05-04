npcdens_demo_source_utils <- function() {
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

npcdens_demo_env <- function(name, default = "") {
  value <- Sys.getenv(name, default)
  if (length(value) != 1L) value <- value[[1L]]
  value
}

npcdens_demo_bool <- function(value, name) {
  value <- trimws(value)
  if (value %in% c("TRUE", "true", "1")) return(TRUE)
  if (value %in% c("FALSE", "false", "0")) return(FALSE)
  stop(name, " must be TRUE/FALSE", call. = FALSE)
}

npcdens_demo_int <- function(value, name, required = TRUE) {
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

npcdens_demo_numeric_field <- function(x) {
  if (is.null(x) || !length(x)) return(NA_real_)
  suppressWarnings(as.numeric(x[[1L]]))
}

npcdens_demo_matrix_path <- function() {
  explicit <- npcdens_demo_env("NP_DEMO_MATRIX")
  if (nzchar(explicit)) return(explicit)
  tier <- npcdens_demo_env("NP_DEMO_TIER", "smoke")
  src <- npcdens_demo_env("NP_DEMO_SRC")
  candidates <- c(
    if (nzchar(src)) file.path(src, "..", "inst", "demo_matrices",
                               paste0("npcdens-", tier, ".csv")) else "",
    system.file("demo_matrices", paste0("npcdens-", tier, ".csv"),
                package = "npRmpi")
  )
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (!length(candidates)) {
    stop("Could not locate npcdens demo matrix for tier=", tier, call. = FALSE)
  }
  candidates[[1L]]
}

npcdens_demo_matrix <- function() {
  path <- npcdens_demo_matrix_path()
  mat <- read.csv(path, stringsAsFactors = FALSE, na.strings = character())
  mat[is.na(mat)] <- ""
  required <- c("family", "case", "tier", "enabled", "default_n", "floor_n",
                "regtype", "bwmethod", "nomad", "degree", "degree.max",
                "bwtype", "why")
  missing <- setdiff(required, names(mat))
  if (length(missing)) {
    stop("npcdens demo matrix missing columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  mat <- mat[mat$enabled == "TRUE", , drop = FALSE]
  cases <- npcdens_demo_env("NP_DEMO_CASES")
  if (!nzchar(cases)) cases <- npcdens_demo_env("NP_DEMO_CASE")
  if (nzchar(cases)) {
    keep <- strsplit(cases, "[[:space:]]+", perl = TRUE)[[1L]]
    keep <- keep[nzchar(keep)]
    mat <- mat[mat$case %in% keep, , drop = FALSE]
  }
  if (!nrow(mat)) stop("npcdens demo matrix selected no enabled rows", call. = FALSE)
  mat
}

npcdens_demo_validate_row <- function(row) {
  if (!identical(row$family, "npcdens")) {
    stop("npcdens family script received family=", row$family, call. = FALSE)
  }
  if (!(row$regtype %in% c("lc", "ll", "lp"))) {
    stop("unsupported npcdens regtype=", row$regtype, call. = FALSE)
  }
  if (!(row$bwmethod %in% c("cv.ls", "cv.ml"))) {
    stop("unsupported npcdens bwmethod=", row$bwmethod, call. = FALSE)
  }
  if (!identical(row$bwtype, "fixed")) {
    stop("npcdens demo matrix currently supports bwtype=fixed only", call. = FALSE)
  }
  nomad <- npcdens_demo_bool(row$nomad, "nomad")
  if (nomad && !identical(row$regtype, "lp")) {
    stop("npcdens nomad rows must use regtype=lp", call. = FALSE)
  }
  if (nomad && nzchar(row$degree)) {
    stop("npcdens nomad rows must leave degree empty for degree search",
         call. = FALSE)
  }
  if (nomad && !nzchar(row$degree.max)) {
    stop("npcdens nomad rows must provide degree.max", call. = FALSE)
  }
  if (!nomad && identical(row$regtype, "lp") && !nzchar(row$degree)) {
    stop("npcdens lp rows must provide degree when nomad=FALSE", call. = FALSE)
  }
  if (!nomad && !identical(row$regtype, "lp") && nzchar(row$degree)) {
    stop("npcdens degree must be empty unless regtype=lp and nomad=FALSE",
         call. = FALSE)
  }
  invisible(TRUE)
}

npcdens_demo_run_row <- function(row, mode) {
  npcdens_demo_validate_row(row)
  default_n <- npcdens_demo_int(row$default_n, "default_n")
  floor_n <- npcdens_demo_int(row$floor_n, "floor_n")
  nomad <- npcdens_demo_bool(row$nomad, "nomad")
  degree <- if (identical(row$regtype, "lp") && !nomad) {
    npcdens_demo_int(row$degree, "degree")
  } else {
    NA_integer_
  }
  degree.max <- if (nomad) npcdens_demo_int(row$degree.max, "degree.max") else NA_integer_
  n <- np_demo_n(default_n, floor = floor_n)

  set.seed(42)
  rho <- 0.25
  data <- MASS::mvrnorm(n = n, mu = c(0, 0), Sigma = matrix(c(1, rho, rho, 1), 2, 2))
  mydat <- data.frame(x = data[, 2], y = data[, 1])
  rm(data)

  bw_args <- list(formula = y ~ x, data = mydat, regtype = row$regtype,
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
    t <- system.time(mpi.bcast.cmd(bw <- do.call(npcdensbw, bw_args),
                                    caller.execute = TRUE))
    t <- t + system.time(mpi.bcast.cmd(model <- npcdens(bws = bw),
                                       caller.execute = TRUE))
  } else {
    t <- system.time(bw <- do.call(npcdensbw, bw_args))
    t <- t + system.time(model <- npcdens(bws = bw))
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
                 selected.degree = npcdens_demo_numeric_field(bw$degree),
                 bwtype = row$bwtype,
                 why = row$why,
                 fitted.length = length(model$condens),
                 objective = npcdens_demo_numeric_field(bw$fval),
                 num.feval = npcdens_demo_numeric_field(bw$num.feval),
                 num.feval.fast = npcdens_demo_numeric_field(bw$num.feval.fast),
                 nomad.time = npcdens_demo_numeric_field(bw$nomad.time),
                 powell.time = npcdens_demo_numeric_field(bw$powell.time))
}

npcdens_demo_run_matrix <- function(mode) {
  npcdens_demo_source_utils()
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("MASS is required for the npcdens demo", call. = FALSE)
  }
  mat <- npcdens_demo_matrix()
  for (i in seq_len(nrow(mat))) {
    row <- mat[i, , drop = FALSE]
    cat("[npcdens][", mode, "] ", row$case, "\n", sep = "")
    npcdens_demo_run_row(row, mode)
  }
  invisible(TRUE)
}
