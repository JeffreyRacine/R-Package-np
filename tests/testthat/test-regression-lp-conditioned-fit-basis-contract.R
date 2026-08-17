library(npRmpi)

.conditioned_fit_test_env <- new.env(parent = emptyenv())

.ensure_conditioned_fit_pool <- function() {
  if (!isTRUE(.conditioned_fit_test_env$started)) {
    skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
    .conditioned_fit_test_env$started <- TRUE
    withr::defer({
      if (isTRUE(.conditioned_fit_test_env$started)) {
        try(close_mpi_slaves(force = TRUE), silent = TRUE)
        .conditioned_fit_test_env$started <- FALSE
      }
    }, envir = testthat::teardown_env())
  }
}

locate_conditioned_fit_source <- function(file) {
  candidates <- c(
    test_path("..", "..", "src", file),
    test_path("..", "..", "..", "src", file),
    file.path(Sys.getenv("R_PACKAGE_DIR", ""), "src", file),
    file.path(Sys.getenv("R_PACKAGE_SOURCE", ""), "src", file),
    file.path(getwd(), "src", file),
    file.path(getwd(), "..", "src", file)
  )
  candidates <- unique(candidates[nzchar(candidates)])
  hits <- candidates[file.exists(candidates)]
  if (!length(hits)) NULL else hits[[1L]]
}

test_that("calendar-scale conditional LP fits use stable global coordinates", {
  .ensure_conditioned_fit_pool()
  old <- options(np.messages = FALSE, np.tree = FALSE, np.largeh = FALSE)
  on.exit(options(old), add = TRUE)

  data("Italy", package = "npRmpi")
  xdat <- data.frame(year = as.numeric(as.character(Italy$year)))
  ydat <- data.frame(gdp = Italy$gdp)
  bandwidth <- c(0.57029458300477509, 4.6970422009288111)

  fit <- function(stable) {
    state <- npcdensbw(
      xdat = xdat, ydat = ydat, bws = bandwidth,
      bandwidth.compute = FALSE, bwscaling = FALSE,
      bwmethod = "cv.ml", bwtype = "fixed",
      regtype = "lp", degree = 3L,
      bernstein.basis = stable
    )
    npcdens(bws = state)$condens
  }

  raw <- fit(FALSE)
  stable <- fit(TRUE)
  expect_true(all(is.finite(raw)))
  expect_equal(raw, stable, tolerance = 2e-11)
})

test_that("ordinary and all-large LP fits share the conditioned basis owner", {
  .ensure_conditioned_fit_pool()
  old <- options(np.messages = FALSE, np.tree = FALSE, np.largeh = FALSE)
  on.exit(options(old), add = TRUE)

  set.seed(20260817L)
  n <- 64L
  xdat <- data.frame(x = 20 + seq(0, 1, length.out = n))
  centered <- xdat$x - 20.5
  ydat <- 0.5 + 0.2 * centered - 0.3 * centered^2 +
    0.1 * centered^3 + rnorm(n, sd = 0.03)

  fit <- function(largeh, stable) {
    options(np.largeh = largeh)
    state <- npregbw(
      xdat = xdat, ydat = ydat, bws = 100,
      bandwidth.compute = FALSE, bwscaling = FALSE,
      bwmethod = "cv.ls", bwtype = "fixed",
      regtype = "lp", degree = 3L,
      bernstein.basis = stable,
      ckertype = "epanechnikov"
    )
    npreg(bws = state, exdat = xdat, gradients = TRUE, se = TRUE)
  }

  for (largeh in c(FALSE, TRUE)) {
    raw <- fit(largeh, FALSE)
    stable <- fit(largeh, TRUE)

    expect_true(all(is.finite(c(raw$mean, raw$grad, raw$merr, raw$gerr))))
    expect_equal(raw$mean, stable$mean, tolerance = 2e-12,
                 info = paste("mean largeh", largeh))
    expect_equal(raw$grad, stable$grad, tolerance = 2e-12,
                 info = paste("gradient largeh", largeh))
    expect_equal(raw$merr, stable$merr, tolerance = 2e-12,
                 info = paste("mean se largeh", largeh))
    expect_equal(raw$gerr, stable$gerr, tolerance = 2e-12,
                 info = paste("gradient se largeh", largeh))
  }
})

test_that("design validation is units neutral but preserves rank errors", {
  calendar <- data.frame(year = 1940 + seq(0, 40, length.out = 80L))
  expect_invisible(npRmpi:::npCheckRegressionDesignCondition(
    reg.code = npRmpi:::REGTYPE_LP,
    xcon = calendar,
    basis = "glp",
    degree = 3L,
    bernstein.basis = FALSE,
    where = "test"
  ))

  collinear <- data.frame(x1 = seq_len(20L), x2 = seq_len(20L))
  expect_error(
    npRmpi:::npCheckRegressionDesignCondition(
      reg.code = npRmpi:::REGTYPE_LP,
      xcon = collinear,
      basis = "glp",
      degree = c(1L, 1L),
      bernstein.basis = FALSE,
      where = "test"
    ),
    "rank deficient",
    fixed = TRUE
  )
})

test_that("fit basis selection is one representation-neutral native seam", {
  jksum_file <- locate_conditioned_fit_source("jksum.c")
  basis_file <- locate_conditioned_fit_source("jksum_lp_basis.c")
  skip_if(is.null(jksum_file), "source file src/jksum.c unavailable")
  skip_if(is.null(basis_file),
          "source file src/jksum_lp_basis.c unavailable")

  jksum <- paste(readLines(jksum_file, warn = FALSE), collapse = "\n")
  basis <- paste(readLines(basis_file, warn = FALSE), collapse = "\n")
  helper_start <- regexpr(
    "static int np_glp_fit_basis_prepare(", jksum, fixed = TRUE
  )[[1L]]
  helper_stop <- regexpr(
    "typedef struct {\n  int ncon;", jksum, fixed = TRUE
  )[[1L]]
  general_start <- regexpr(
    "static SEXP np_regression_general_lp_fit_execute(void *data)",
    jksum, fixed = TRUE
  )[[1L]]
  general_stop <- regexpr(
    "static int np_regression_general_lp_fit(", jksum, fixed = TRUE
  )[[1L]]

  expect_gt(helper_start, 0L)
  expect_gt(helper_stop, helper_start)
  expect_gt(general_start, 0L)
  expect_gt(general_stop, general_start)

  helper <- substr(jksum, helper_start, helper_stop - 1L)
  general <- substr(jksum, general_start, general_stop - 1L)
  all_large <- substr(jksum, general_stop, nchar(jksum))

  expect_match(helper, "np_lp_basis_requires_conditioning(", fixed = TRUE)
  expect_match(helper, "NP_GLP_RAW_BASIS_MIN_RCOND", fixed = TRUE)
  expect_match(helper, "np_glp_basis_ctx_array_prepare(", fixed = TRUE)
  expect_match(helper, "np_glp_fill_basis_raw_train(", fixed = TRUE)
  expect_match(helper, "np_glp_fill_basis_train(", fixed = TRUE)
  expect_false(grepl("vector_glp_degree_extern", helper, fixed = TRUE))
  expect_match(general, "np_glp_fit_basis_prepare(", fixed = TRUE)
  expect_match(all_large, "np_glp_fit_basis_prepare_execute", fixed = TRUE)
  expect_false(grepl("NPRegressionAllLargeBernsteinBasisCall", jksum,
                     fixed = TRUE))
  expect_false(grepl("np_regression_alllarge_bernstein_basis_execute", jksum,
                     fixed = TRUE))

  expect_match(basis, "if(contiguous)", fixed = TRUE)
  expect_match(basis, "F77_CALL(dsyrk)", fixed = TRUE)
  expect_match(basis, "source[i][row]*source[j][row]", fixed = TRUE)
  expect_false(grepl("n*n", helper, fixed = TRUE))
})
