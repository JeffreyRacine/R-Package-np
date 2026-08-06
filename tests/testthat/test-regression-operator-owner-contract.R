locate_regression_operator_source <- function() {
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

regression_operator_owner_region <- function(source) {
  start <- grep(
    "^int kernel_estimate_regression_categorical_tree_np\\(", source
  )
  stop <- grep("^static int np_conditional_indicator_row_core\\(", source)
  stopifnot(length(start) == 1L, length(stop) == 1L, stop > start)
  paste(source[seq.int(start, stop - 1L)], collapse = "\n")
}

test_that("regression fit has one post-validation operator owner", {
  path <- locate_regression_operator_source()
  skip_if(is.null(path), "package sources unavailable")
  source <- readLines(path, warn = FALSE)
  region <- regression_operator_owner_region(source)

  allocations <- gregexpr(
    "operator = num_reg_total > 0 ? (int *)malloc(bytes) : NULL;",
    region,
    fixed = TRUE
  )[[1L]]
  expect_equal(sum(allocations > 0L), 1L)
  allocation <- allocations[allocations > 0L][[1L]]
  scalar_return <- regexpr(
    "lp_engine_est == NP_LP_ENGINE_SCALAR", region, fixed = TRUE
  )
  expect_gt(scalar_return, 0L)
  expect_gt(allocation, scalar_return)
  expect_match(region, "fit_owner.operator = operator;", fixed = TRUE)
  expect_false(grepl("free(operator);", region, fixed = TRUE))
  expect_equal(
    sum(gregexpr(
      "np_regression_fit_owner_clear(&fit_owner);",
      region,
      fixed = TRUE
    )[[1L]] > 0L),
    1L
  )
  owner_clear_start <- grep(
    "^static void np_regression_fit_owner_clear\\(", source
  )
  owner_clear_stop <- grep(
    "^static void np_regression_fit_owner_cleanup\\(", source
  )
  expect_length(owner_clear_start, 1L)
  expect_length(owner_clear_stop, 1L)
  owner_clear <- paste(
    source[owner_clear_start:(owner_clear_stop - 1L)], collapse = "\n"
  )
  expect_equal(
    sum(gregexpr("free(owner->operator);", owner_clear, fixed = TRUE)[[1L]] > 0L),
    1L
  )
})

test_that("scalar and general beta regression repeats are exact", {
  set.seed(20260804)
  training <- data.frame(
    x = sort(runif(48L, 0.04, 0.96)),
    z = runif(48L, 0.05, 0.95)
  )
  response <- sin(2 * pi * training$x) + 0.35 * training$z
  evaluation <- training[c(3L, 11L, 22L, 37L, 46L), , drop = FALSE]

  fit_once <- function(degree) {
    fit <- npreg(
      bws = c(0.21, 0.24),
      txdat = training,
      tydat = response,
      exdat = evaluation,
      gradients = TRUE,
      regtype = "lp",
      degree = rep.int(degree, 2L),
      bernstein.basis = FALSE,
      ckertype = "beta",
      ckerorder = 6L,
      ckerbound = "fixed",
      ckerlb = c(0, 0),
      ckerub = c(1, 1)
    )
    list(
      mean = fitted(fit),
      se = se(fit),
      gradients = gradients(fit),
      gradient_se = fit$gerr
    )
  }

  for (degree in c(0L, 1L)) {
    reference <- fit_once(degree)
    for (iteration in seq_len(12L))
      expect_identical(fit_once(degree), reference)
  }
})
