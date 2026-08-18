.conditional_lp_conditioning_test_env <- new.env(parent = emptyenv())

.ensure_conditional_lp_conditioning_pool <- function() {
  if (!isTRUE(.conditional_lp_conditioning_test_env$started)) {
    skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
    .conditional_lp_conditioning_test_env$started <- TRUE
    withr::defer({
      if (isTRUE(.conditional_lp_conditioning_test_env$started)) {
        try(close_mpi_slaves(force = TRUE), silent = TRUE)
        .conditional_lp_conditioning_test_env$started <- FALSE
      }
    }, envir = testthat::teardown_env())
  }
}

locate_conditional_lp_source <- function(file) {
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

test_that("raw conditional LP searches are invariant to polynomial coordinates", {
  .ensure_conditional_lp_conditioning_pool()
  old <- options(np.messages = FALSE, np.tree = FALSE, np.largeh = FALSE)
  on.exit(options(old), add = TRUE)

  data("Italy", package = "npRmpi")
  italy <- within(Italy, year <- as.numeric(as.character(year)))

  density <- npcdens(
    gdp ~ year, data = italy,
    regtype = "lp", degree = 1L, nmulti = 1L
  )
  distribution <- npcdist(
    gdp ~ year, data = italy,
    regtype = "lp", degree = 1L, nmulti = 1L
  )

  expect_true(is.finite(density$bws$fval))
  expect_true(is.finite(distribution$bws$fval))
  expect_equal(density$bws$fval, -2622.1740754235052, tolerance = 1e-10)
  expect_equal(distribution$bws$fval, 0.088110070138586946,
               tolerance = 1e-12)
})

test_that("conditional influence objectives are basis-family neutral", {
  .ensure_conditional_lp_conditioning_pool()
  old <- options(np.messages = FALSE, np.tree = FALSE, np.largeh = FALSE)
  on.exit(options(old), add = TRUE)

  data("Italy", package = "npRmpi")
  xdat <- data.frame(year = as.numeric(as.character(Italy$year)))
  ydat <- data.frame(gdp = Italy$gdp)
  density_bandwidths <- list(
    c(2.1291913099237134, 0.54853849744276473),
    c(4.2798952993925949, 0.60863958901345139),
    c(4.6970422009288111, 0.57029458300477509),
    c(6.1492474834365085, 0.64670985511180268),
    c(7.5557206581087621, 0.59213538803902854)
  )
  distribution_bandwidths <- list(
    c(3.2974184460991904, 0.36695247958275035),
    c(5.3363521464555861, 0.52753928548220452),
    c(6.018324870130769, 0.55292169796488),
    c(7.6327016627794873, 0.62504080125894779),
    c(10.212704973037431, 0.69714228399430267)
  )

  evaluate <- function(constructor, evaluator, method, bandwidth, degree,
                       bernstein) {
    state <- constructor(
      xdat = xdat, ydat = ydat, bws = bandwidth,
      bandwidth.compute = FALSE, bwscaling = FALSE,
      bwmethod = method, bwtype = "fixed",
      regtype = "lp", degree = degree,
      bernstein.basis = bernstein
    )
    as.numeric(evaluator(xdat = xdat, ydat = ydat, bws = state)$objective)
  }

  for (degree in 1:5) {
    density_raw <- evaluate(
      npcdensbw, npRmpi:::.npcdensbw_eval_only, "cv.ml",
      density_bandwidths[[degree]], degree, FALSE
    )
    density_graded <- evaluate(
      npcdensbw, npRmpi:::.npcdensbw_eval_only, "cv.ml",
      density_bandwidths[[degree]], degree, TRUE
    )
    distribution_raw <- evaluate(
      npcdistbw, npRmpi:::.npcdistbw_eval_only, "cv.ls",
      distribution_bandwidths[[degree]], degree, FALSE
    )
    distribution_graded <- evaluate(
      npcdistbw, npRmpi:::.npcdistbw_eval_only, "cv.ls",
      distribution_bandwidths[[degree]], degree, TRUE
    )

    expect_equal(density_raw, density_graded, tolerance = 5e-11,
                 info = paste("density", degree))
    expect_equal(distribution_raw, distribution_graded, tolerance = 5e-11,
                 info = paste("distribution", degree))
  }
})

test_that("conditional LP conditioning is translation stable and topology general", {
  .ensure_conditional_lp_conditioning_pool()
  old <- options(np.messages = FALSE, np.tree = FALSE, np.largeh = FALSE)
  on.exit(options(old), add = TRUE)

  set.seed(2026080201L)
  n <- 96L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  translated <- transform(x, x1 = x1 + 1950, x2 = x2 - 720)
  y <- data.frame(
    y = 0.3 + sin(2 * pi * x$x1) * cos(pi * x$x2) +
      rnorm(n, sd = 0.08)
  )

  objective <- function(xdat, basis, bernstein) {
    state <- npcdensbw(
      xdat = xdat, ydat = y, bws = c(0.28, 0.31, 0.27),
      bandwidth.compute = FALSE, bwscaling = FALSE,
      bwmethod = "cv.ls", bwtype = "fixed",
      regtype = "lp", basis = basis, degree = c(2L, 2L),
      bernstein.basis = bernstein,
      cxkertype = "gaussian", cxkerorder = 2L,
      cykertype = "gaussian", cykerorder = 2L
    )
    as.numeric(npRmpi:::.npcdensbw_eval_only(xdat, y, state)$objective)
  }

  for (basis in c("glp", "additive", "tensor")) {
    raw <- objective(x, basis, FALSE)
    graded <- objective(x, basis, TRUE)
    raw_translated <- objective(translated, basis, FALSE)
    graded_translated <- objective(translated, basis, TRUE)

    expect_equal(raw, graded, tolerance = 2e-12, info = basis)
    expect_equal(raw_translated, graded_translated,
                 tolerance = 2e-12, info = basis)
    expect_equal(raw_translated, raw, tolerance = 2e-10, info = basis)
  }
})

test_that("conditional LP influence coordinates are one cached global policy", {
  jksum_file <- locate_conditional_lp_source("jksum.c")
  basis_file <- locate_conditional_lp_source("jksum_lp_basis.c")
  skip_if(is.null(jksum_file), "source file src/jksum.c unavailable")
  skip_if(is.null(basis_file),
          "source file src/jksum_lp_basis.c unavailable")

  jksum <- paste(readLines(jksum_file, warn = FALSE), collapse = "\n")
  basis <- paste(readLines(basis_file, warn = FALSE), collapse = "\n")
  start <- regexpr(
    "static int np_glp_cv_cache_prepare(const int lp_engine,",
    jksum, fixed = TRUE
  )[[1L]]
  stop <- regexpr(
    "static int np_glp_cv_cache_prepare_original_order(const int *ipt)",
    jksum, fixed = TRUE
  )[[1L]]
  expect_gt(start, 0L)
  expect_gt(stop, start)
  prepare <- substr(jksum, start, stop - 1L)

  expect_match(prepare, "np_glp_fit_basis_prepare(", fixed = TRUE)
  expect_match(prepare, "source_basis_stable = use_stable_basis;",
               fixed = TRUE)
  expect_match(prepare, "np_glp_cv_cache.basis = np_glp_cv_cache.source_basis;",
               fixed = TRUE)
  expect_match(prepare, "influence_basis_ready = 1;", fixed = TRUE)
  expect_match(
    prepare,
    "influence_basis_conditioned =\n    np_glp_cv_cache.source_basis_stable;",
    fixed = TRUE
  )
  expect_false(grepl("np_lp_conditioned_basis_prepare(", prepare,
                     fixed = TRUE))
  expect_false(grepl("conditioned_basis", prepare, fixed = TRUE))
  expect_false(grepl("dgeqp3", jksum, fixed = TRUE))
  expect_false(grepl("dorgqr", jksum, fixed = TRUE))
  expect_false(grepl("kernel_c", prepare, fixed = TRUE))
  expect_false(grepl("NPGLPQRDropWorkspace", jksum, fixed = TRUE))
  expect_false(grepl("num_obs*num_obs", prepare, fixed = TRUE))
  expect_false(grepl("bernstein", basis, ignore.case = TRUE))
  expect_false(grepl("kernel", basis, ignore.case = TRUE))
  expect_match(basis, "F77_CALL(dsyrk)", fixed = TRUE)
  expect_false(exists(
    "npWithLocalLinearRawBasisSearchError",
    envir = asNamespace("npRmpi"), inherits = FALSE
  ))

  row_start <- regexpr(
    "static int np_conditional_x_weight_row_stream_core_impl(",
    jksum, fixed = TRUE
  )[[1L]]
  row_stop <- regexpr(
    "static int np_conditional_x_weight_row_stream_core(",
    jksum, fixed = TRUE
  )[[1L]]
  expect_gt(row_start, 0L)
  expect_gt(row_stop, row_start)
  row_owner <- substr(jksum, row_start, row_stop - 1L)
  expect_match(row_owner, "if(!np_glp_cv_cache.ready", fixed = TRUE)
  expect_false(grepl("np_glp_cv_clear_extern();", row_owner, fixed = TRUE))

  cvls_start <- regexpr(
    "static int np_conditional_density_cvls_lp_row_stream(",
    jksum, fixed = TRUE
  )[[1L]]
  cvls_stop <- regexpr(
    "#define NP_CDENS_ADAP_WIDTH3_NOT_BENEFICIAL",
    jksum, fixed = TRUE
  )[[1L]]
  expect_gt(cvls_start, 0L)
  expect_gt(cvls_stop, cvls_start)
  cvls_owner <- substr(jksum, cvls_start, cvls_stop - 1L)
  expect_match(cvls_owner, "np_glp_cv_clear_extern();", fixed = TRUE)

  roots <- unique(c(
    test_path("..", ".."),
    test_path("..", "..", ".."),
    Sys.getenv("R_PACKAGE_DIR", ""),
    Sys.getenv("R_PACKAGE_SOURCE", ""),
    getwd(),
    file.path(getwd(), "..")
  ))
  roots <- roots[nzchar(roots)]
  roots <- roots[file.exists(file.path(roots, "R", "np.condensity.R"))]
  skip_if(!length(roots), "package R sources unavailable")
  conditional_dispatch <- paste(
    readLines(file.path(roots[[1L]], "R", "np.condensity.R"), warn = FALSE),
    readLines(file.path(roots[[1L]], "R", "np.condistribution.R"),
              warn = FALSE),
    readLines(file.path(roots[[1L]], "R", "util.R"), warn = FALSE),
    collapse = "\n"
  )
  expect_false(grepl("keep_local_raw_degree1_cvls", conditional_dispatch,
                     fixed = TRUE))
  expect_false(grepl("npIsRawDegreeOneConditional", conditional_dispatch,
                     fixed = TRUE))
})
