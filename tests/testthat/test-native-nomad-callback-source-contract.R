np_extract_c_function_body <- function(source_text, name) {
  pattern <- paste0("(?m)^\\s*static\\s+[^\\n;{]*\\b", name, "\\s*\\(")
  locs <- gregexpr(pattern, source_text, perl = TRUE)[[1L]]
  if (identical(locs, -1L))
    stop("missing C function ", name, call. = FALSE)

  chars <- strsplit(source_text, "", fixed = TRUE)[[1L]]
  definition_locs <- integer()
  for (loc in locs) {
    open_try <- regexpr("\\{", substring(source_text, loc), perl = TRUE)[[1L]]
    semi_try <- regexpr(";", substring(source_text, loc), perl = TRUE)[[1L]]
    if (open_try > 0L && (semi_try <= 0L || open_try < semi_try))
      definition_locs <- c(definition_locs, loc)
  }
  if (!length(definition_locs))
    stop("missing C function definition ", name, call. = FALSE)
  if (length(definition_locs) > 1L)
    stop("multiple C function definitions ", name, call. = FALSE)

  loc <- definition_locs[[1L]]
  open <- NA_integer_
  for (i in seq.int(loc, length(chars))) {
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
  source_file <- file.path(np_namespace_hygiene_root(), "src", "np.c")
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
    "np_density_conditional_nomad_shadow_eval_native_raw",
    "np_regression_native_decode_eval_bw",
    "np_density_prepared_context_eval",
    "np_distribution_prepared_context_eval",
    "np_distribution_conditional_nomad_native_eval_once",
    "np_udens_native_decode_eval_bw",
    "np_udist_native_decode_eval_bw",
    "np_cdist_native_decode_eval_bw",
    "np_cdens_native_search_callback",
    "np_regression_native_search_callback",
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
    body <- np_extract_c_function_body(source_text, fun)
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
  expect_false(grepl(
    "\\bnp_density_nomad_native_eval_once\\b",
    source_text,
    perl = TRUE
  ))
  expect_false(grepl(
    "\\bnp_regression_nomad_native_eval_once\\b",
    source_text,
    perl = TRUE
  ))
  expect_false(grepl(
    "\\bnp_distribution_nomad_native_eval_once\\b",
    source_text,
    perl = TRUE
  ))
  expect_false(grepl(
    "NP_NOMAD_CALLBACK_CALLOC",
    np_extract_c_function_body(source_text, "np_udens_native_search_callback"),
    fixed = TRUE
  ))
  init_text <- paste(readLines(
    file.path(np_namespace_hygiene_root(), "src", "np_init.c"),
    warn = FALSE
  ), collapse = "\n")
  expect_false(grepl(
    "C_np_density_nomad_native_fixed_eval",
    init_text,
    fixed = TRUE
  ))
})

test_that("native NOMAD categorical coordinates use one canonical decoder", {
  source_file <- file.path(np_namespace_hygiene_root(), "src", "np.c")
  expect_true(file.exists(source_file))
  source_text <- paste(readLines(source_file, warn = FALSE), collapse = "\n")
  helper <- np_extract_c_function_body(
    source_text, "np_nomad_decode_categorical_bandwidth"
  )

  expect_match(helper, "raw_point/1.0e4", fixed = TRUE)
  expect_match(helper, "external/ncatfac", fixed = TRUE)
  for (fun in c(
    "np_regression_native_decode_eval_bw",
    "np_udens_native_decode_eval_bw",
    "np_udist_native_decode_eval_bw",
    "np_cdist_native_decode_eval_bw",
    "np_cdens_native_search_callback"
  )) {
    expect_match(
      np_extract_c_function_body(source_text, fun),
      "np_nomad_decode_categorical_bandwidth(",
      fixed = TRUE,
      info = fun
    )
  }
})
