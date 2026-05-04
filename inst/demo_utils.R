np_demo_n <- function(default,
                      floor = 20L,
                      exact_env = "NP_DEMO_N",
                      frac_env = "NP_DEMO_N_FRAC") {
  validate_scalar_positive <- function(x, name, integer = FALSE) {
    if (length(x) != 1L || is.na(x) || !is.finite(x) || x <= 0) {
      stop(name, " must be a positive ", if (integer) "integer" else "number", call. = FALSE)
    }
    if (integer && x != as.integer(x)) {
      stop(name, " must be a positive integer", call. = FALSE)
    }
    invisible(TRUE)
  }

  validate_scalar_positive(default, "default", integer = TRUE)
  validate_scalar_positive(floor, "floor", integer = TRUE)

  parse_exact <- function(value, env_name) {
    parsed <- suppressWarnings(as.numeric(value))
    if (length(parsed) != 1L || is.na(parsed) || !is.finite(parsed) ||
        parsed <= 0 || parsed != as.integer(parsed)) {
      stop(env_name, " must be a positive integer", call. = FALSE)
    }
    as.integer(parsed)
  }

  parse_fraction <- function(value, env_name) {
    parsed <- suppressWarnings(as.numeric(value))
    if (length(parsed) != 1L || is.na(parsed) || !is.finite(parsed) || parsed <= 0) {
      stop(env_name, " must be a positive number", call. = FALSE)
    }
    parsed
  }

  exact <- Sys.getenv(exact_env, "")
  if (nzchar(exact)) {
    return(as.integer(max(as.integer(floor), parse_exact(exact, exact_env))))
  }

  fraction <- Sys.getenv(frac_env, "")
  if (nzchar(fraction)) {
    n <- ceiling(as.integer(default) * parse_fraction(fraction, frac_env))
    return(as.integer(max(as.integer(floor), n)))
  }

  as.integer(default)
}

np_demo_result <- function(demo, mode, n, default_n, elapsed, ...) {
  fields <- list(
    demo = demo,
    mode = mode,
    n = as.integer(n),
    default_n = as.integer(default_n),
    elapsed = sprintf("%.6f", as.numeric(elapsed))
  )
  fields <- c(fields, list(...))

  encode <- function(x) {
    if (length(x) == 0L || is.null(x)) return("NA")
    x <- x[[1L]]
    if (is.na(x)) return("NA")
    if (is.numeric(x)) return(sprintf("%.12g", x))
    gsub("[[:space:]]+", "_", as.character(x), perl = TRUE)
  }

  cat(
    "DEMO_RESULT ",
    paste(sprintf("%s=%s", names(fields), vapply(fields, encode, character(1L))), collapse = " "),
    "\n",
    sep = ""
  )
}
