locate_categorical_profile_owner_source <- function() {
  roots <- unique(c(
    testthat::test_path("..", ".."),
    testthat::test_path("..", "..", ".."),
    Sys.getenv("R_PACKAGE_DIR", ""),
    Sys.getenv("R_PACKAGE_SOURCE", ""),
    getwd(),
    file.path(getwd(), "..")
  ))
  roots <- roots[nzchar(roots)]
  paths <- file.path(roots, "src", "jksum.c")
  paths <- paths[file.exists(paths)]
  if (!length(paths))
    return(NULL)
  paths[[1L]]
}

categorical_profile_body <- function(source, helper) {
  start_pattern <- paste0("^static int np_", helper, "_body\\(")
  finish_pattern <- paste0("^static SEXP np_", helper, "_execute\\(")
  start <- grep(start_pattern, source)
  finish <- grep(finish_pattern, source)
  stopifnot(length(start) == 1L, length(finish) == 1L, finish > start)
  paste(source[seq.int(start, finish - 1L)], collapse = "\n")
}

fixed_count <- function(text, pattern) {
  locations <- gregexpr(pattern, text, fixed = TRUE)[[1L]]
  sum(locations > 0L)
}

test_that("five categorical-profile helpers share one complete owner", {
  path <- locate_categorical_profile_owner_source()
  skip_if(is.null(path), "package sources unavailable")
  source <- readLines(path, warn = FALSE)
  all_source <- paste(source, collapse = "\n")
  compact <- gsub("[[:space:]]+", " ", all_source)

  expect_match(
    compact,
    paste0(
      "typedef struct { int *indices[NP_CATEGORICAL_PROFILE_OWNER_CAPACITY]; ",
      "double **matrices[NP_CATEGORICAL_PROFILE_OWNER_CAPACITY]; ",
      "double *vectors[NP_CATEGORICAL_PROFILE_OWNER_CAPACITY];"
    ),
    fixed = TRUE
  )
  expect_match(
    compact,
    "NP_CATEGORICAL_PROFILE_OWNER_CAPACITY = 8",
    fixed = TRUE
  )
  expect_match(
    compact,
    paste0(
      "if(jump) np_categorical_profile_owner_clear(",
      "(NPCategoricalProfileOwner *)data);"
    ),
    fixed = TRUE
  )

  helpers <- c(
    "regression_categorical_profile_cv",
    "regression_categorical_profile_fit",
    "density_categorical_profile_fit",
    "density_categorical_profile_cv",
    "conditional_categorical_profile_fit"
  )
  expected <- list(
    regression_categorical_profile_cv = c(index = 2L, matrix = 2L,
                                          vector = 4L, kernel = 1L),
    regression_categorical_profile_fit = c(index = 4L, matrix = 4L,
                                           vector = 6L, kernel = 1L),
    density_categorical_profile_fit = c(index = 4L, matrix = 4L,
                                        vector = 2L, kernel = 1L),
    density_categorical_profile_cv = c(index = 3L, matrix = 2L,
                                       vector = 3L, kernel = 2L),
    conditional_categorical_profile_fit = c(index = 8L, matrix = 8L,
                                            vector = 4L, kernel = 2L)
  )

  for (helper in helpers) {
    body <- categorical_profile_body(source, helper)
    counts <- expected[[helper]]
    expect_equal(
      fixed_count(body, "np_categorical_profile_owner_take_index("),
      unname(counts[["index"]]),
      info = helper
    )
    expect_equal(
      fixed_count(body, "np_categorical_profile_owner_take_matrix("),
      unname(counts[["matrix"]]),
      info = helper
    )
    expect_equal(
      fixed_count(body, "np_categorical_profile_owner_take_vector("),
      unname(counts[["vector"]]),
      info = helper
    )
    expect_equal(
      fixed_count(body, "kernel_weighted_sum_np("),
      unname(counts[["kernel"]]),
      info = helper
    )
    expect_equal(
      fixed_count(body, "np_categorical_profile_owner_clear(owner);"),
      1L,
      info = helper
    )
    expect_false(grepl("if\\([^\\n]*!= NULL\\) free", body),
                 info = helper)
    expect_false(grepl("num_obs(_train|_eval)?\\s*\\*\\s*num_obs",
                       body), info = helper)

    protected_call <- paste0(
      "R_UnwindProtect( np_", helper, "_execute, &execution, ",
      "np_categorical_profile_owner_cleanup, &execution.owner, NULL);"
    )
    expect_match(compact, protected_call, fixed = TRUE, info = helper)
  }
})

test_that("categorical-profile owner capacity failures cannot fall back", {
  path <- locate_categorical_profile_owner_source()
  skip_if(is.null(path), "package sources unavailable")
  compact <- gsub(
    "[[:space:]]+", " ",
    paste(readLines(path, warn = FALSE), collapse = "\n")
  )

  for (resource in c("index", "matrix", "vector")) {
    expect_match(
      compact,
      paste0(
        'error("internal categorical-profile ', resource,
        '-owner capacity exceeded")'
      ),
      fixed = TRUE
    )
  }
  expect_match(
    compact,
    paste0(
      "if(owner->index_count >= NP_CATEGORICAL_PROFILE_OWNER_CAPACITY) { ",
      "free(index); error("
    ),
    fixed = TRUE
  )
  expect_match(
    compact,
    paste0(
      "if(owner->matrix_count >= NP_CATEGORICAL_PROFILE_OWNER_CAPACITY) { ",
      "free_tmat(matrix); error("
    ),
    fixed = TRUE
  )
  expect_match(
    compact,
    paste0(
      "if(owner->vector_count >= NP_CATEGORICAL_PROFILE_OWNER_CAPACITY) { ",
      "free(vector); error("
    ),
    fixed = TRUE
  )
})

