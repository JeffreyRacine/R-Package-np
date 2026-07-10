npRmpi_extract_c_function_body <- function(source_text, name) {
  loc <- regexpr(paste0("\\b", name, "\\s*\\("), source_text, perl = TRUE)
  if (loc[[1L]] <= 0L)
    stop("missing C function ", name, call. = FALSE)

  chars <- strsplit(source_text, "", fixed = TRUE)[[1L]]
  open <- NA_integer_
  for (i in seq.int(loc[[1L]], length(chars))) {
    if (identical(chars[[i]], "{")) {
      open <- i
      break
    }
  }
  if (is.na(open))
    stop("missing C body for ", name, call. = FALSE)

  depth <- 0L
  close <- NA_integer_
  for (i in seq.int(open, length(chars))) {
    if (identical(chars[[i]], "{")) {
      depth <- depth + 1L
    } else if (identical(chars[[i]], "}")) {
      depth <- depth - 1L
      if (identical(depth, 0L)) {
        close <- i
        break
      }
    }
  }
  if (is.na(close))
    stop("unterminated C body for ", name, call. = FALSE)

  paste0(chars[open:close], collapse = "")
}

test_that("native NOMAD C callback path does not call R API or longjmp helpers", {
  source_file <- file.path(npRmpi_namespace_hygiene_root(), "src", "np.c")
  expect_true(file.exists(source_file))
  source_text <- paste(readLines(source_file, warn = FALSE), collapse = "\n")

  expect_match(source_text, "return\\s+calloc\\s*\\(", perl = TRUE)
  expect_false(grepl(
    "#define\\s+NP_NOMAD_CALLBACK_CALLOC[^\\n]*R_Calloc",
    source_text,
    perl = TRUE
  ))

  callback_path <- c(
    "bwmfunc_wrapper",
    "np_regression_shadow_native_decode_bw",
    "np_density_conditional_nomad_shadow_eval_native_raw",
    "np_density_nomad_native_eval_once",
    "np_distribution_nomad_native_eval_once",
    "np_distribution_conditional_nomad_native_eval_once",
    "np_udens_native_decode_eval_bw",
    "np_udist_native_decode_eval_bw",
    "np_cdist_native_decode_eval_bw",
    "np_regression_shadow_native_search_callback",
    "np_cdens_native_search_callback",
    "np_udens_native_search_callback",
    "np_udist_native_search_callback",
    "np_cdist_native_search_callback"
  )

  forbidden <- c(
    "\\bR_Calloc\\s*\\(",
    "\\bR_Realloc\\s*\\(",
    "\\bR_alloc\\s*\\(",
    "\\bR_Free\\s*\\(",
    "\\bRf_[A-Za-z0-9_]+\\s*\\(",
    "\\berror\\s*\\(",
    "\\bwarning\\s*\\(",
    "\\bPROTECT\\s*\\(",
    "\\bUNPROTECT\\s*\\(",
    "\\ballocVector\\s*\\(",
    "\\bScalar[A-Za-z0-9_]*\\s*\\(",
    "\\bmkString\\s*\\(",
    "\\bR_GetCCallable\\s*\\(",
    "\\bGetOption\\s*\\(",
    "\\bRprintf\\s*\\(",
    "\\bREprintf\\s*\\("
  )

  violations <- character()
  for (fun in callback_path) {
    body <- npRmpi_extract_c_function_body(source_text, fun)
    hits <- forbidden[vapply(forbidden, grepl, logical(1L), x = body, perl = TRUE)]
    if (length(hits)) {
      violations <- c(
        violations,
        sprintf("%s: %s", fun, paste(hits, collapse = ", "))
      )
    }
  }

  if (length(violations))
    stop(paste(violations, collapse = "\n"), call. = FALSE)
  expect_length(violations, 0L)
})
