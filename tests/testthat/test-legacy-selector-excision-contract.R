locate_legacy_selector_sources <- function() {
  roots <- c(
    testthat::test_path("..", ".."),
    testthat::test_path("..", "..", ".."),
    Sys.getenv("R_PACKAGE_DIR", ""),
    Sys.getenv("R_PACKAGE_SOURCE", ""),
    getwd(),
    file.path(getwd(), "..")
  )
  roots <- unique(roots[nzchar(roots)])
  roots <- roots[file.exists(file.path(roots, "src", "np.c"))]
  if (!length(roots))
    return(NULL)
  roots[[1L]]
}

test_that("retired density selectors cannot restore a second engine", {
  root <- locate_legacy_selector_sources()
  skip_if(is.null(root), "package sources unavailable")

  headers <- paste(
    readLines(file.path(root, "src", "headers.h"), warn = FALSE),
    collapse = "\n"
  )
  native <- paste(
    readLines(file.path(root, "src", "np.c"), warn = FALSE),
    collapse = "\n"
  )
  registration <- paste(
    readLines(file.path(root, "src", "np_init.c"), warn = FALSE),
    collapse = "\n"
  )
  r_sources <- paste(
    unlist(lapply(
      c(
        "np.density.bw.R",
        "np.density.R",
        "np.condensity.bw.R"
      ),
      function(file) {
        readLines(file.path(root, "R", file), warn = FALSE)
      }
    )),
    collapse = "\n"
  )

  expect_match(headers, "#define BW_RESERVED_LEGACYI 16", fixed = TRUE)
  expect_match(headers, "#define CBW_RESERVED_LEGACYI 22", fixed = TRUE)
  expect_match(headers, "#define DEN_RESERVED_LEGACYI 14", fixed = TRUE)

  retired_tokens <- c(
    "BW_OLDBW",
    "CBW_OLDI",
    "DEN_OLDI",
    "old_bw",
    "old_cdens",
    "old_dens"
  )
  for (token in retired_tokens) {
    expect_false(grepl(token, paste(headers, native), fixed = TRUE))
  }
  expect_false(grepl("old.dens =", r_sources, fixed = TRUE))
  expect_false(grepl("old.cdens =", r_sources, fixed = TRUE))
  expect_false(grepl("compare_old", native, fixed = TRUE))
  expect_false(grepl("do_old", native, fixed = TRUE))
  expect_false(grepl("mkChar(\"old\")", native, fixed = TRUE))
  retired_proof_symbols <- c(
    "C_np_shadow_cv_density_conditional",
    "C_np_shadow_cv_xweights_conditional",
    "C_np_shadow_cv_yrow_conditional",
    "C_np_shadow_cv_xweights_full_conditional",
    "C_np_shadow_cv_distribution_conditional",
    "C_np_shadow_reset_state",
    "np_shadow_"
  )
  for (symbol in retired_proof_symbols) {
    expect_false(grepl(symbol, paste(headers, native, registration), fixed = TRUE))
  }
  expect_match(
    registration,
    "{\"C_np_reset_native_estimator_state\",(DL_FUNC) &C_np_reset_native_estimator_state,0}",
    fixed = TRUE
  )

  expect_match(
    native,
    "C_np_density_bw: reserved legacy selector must be zero",
    fixed = TRUE
  )
  expect_match(
    native,
    "C_np_density_conditional_bw: reserved legacy selector must be zero",
    fixed = TRUE
  )
  expect_match(
    native,
    "C_np_density: reserved legacy selector must be zero",
    fixed = TRUE
  )
  expect_match(
    native,
    "bwmfunc = np_cv_func_density_categorical_ml;",
    fixed = TRUE
  )
  expect_match(
    native,
    "bwmfunc = np_cv_func_density_categorical_ls;",
    fixed = TRUE
  )
  expect_match(
    native,
    "bwmfunc = np_cv_func_con_density_categorical_ml;",
    fixed = TRUE
  )
  expect_match(
    native,
    "bwmfunc = np_cv_func_con_density_categorical_ls_npksum;",
    fixed = TRUE
  )
})
