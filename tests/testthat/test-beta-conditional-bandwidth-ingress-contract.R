locate_beta_conditional_bandwidth_sources <- function() {
  roots <- unique(c(
    Sys.getenv("NP_SOURCE_ROOT", unset = ""),
    normalizePath(file.path(getwd(), "..", ".."), mustWork = FALSE),
    normalizePath(getwd(), mustWork = FALSE)
  ))
  roots <- roots[nzchar(roots)]
  for (root in roots) {
    if (file.exists(file.path(root, "src", "np.c")) &&
        file.exists(file.path(root, "R", "np.condensity.bw.R")) &&
        file.exists(file.path(root, "R", "np.condistribution.bw.R")))
      return(root)
  }
  NULL
}

.beta_conditional_ingress_test_env <- new.env(parent = emptyenv())

.ensure_beta_conditional_ingress_pool <- function() {
  if (!isTRUE(.beta_conditional_ingress_test_env$started)) {
    skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
    .beta_conditional_ingress_test_env$started <- TRUE
    withr::defer({
      if (isTRUE(.beta_conditional_ingress_test_env$started)) {
        try(close_mpi_slaves(force = TRUE), silent = TRUE)
        .beta_conditional_ingress_test_env$started <- FALSE
      }
    }, envir = testthat::teardown_env())
  }
}

test_that("conditional bandwidth ingress marshals one compression contract", {
  root <- locate_beta_conditional_bandwidth_sources()
  skip_if(is.null(root), "package sources unavailable")
  density_r <- paste(
    readLines(file.path(root, "R", "np.condensity.bw.R"), warn = FALSE),
    collapse = "\n"
  )
  distribution_r <- paste(
    readLines(file.path(root, "R", "np.condistribution.bw.R"), warn = FALSE),
    collapse = "\n"
  )
  headers <- paste(
    readLines(file.path(root, "src", "headers.h"), warn = FALSE),
    collapse = "\n"
  )
  ingress <- paste(
    readLines(file.path(root, "src", "np.c"), warn = FALSE),
    collapse = "\n"
  )

  option_literal <- "categorical.compress = npStrictLogicalOption("
  density_hits <- gregexpr(option_literal, density_r, fixed = TRUE)[[1L]]
  distribution_hits <- gregexpr(
    option_literal, distribution_r, fixed = TRUE
  )[[1L]]
  expect_length(density_hits[density_hits > 0L], 3L)
  expect_length(distribution_hits[distribution_hits > 0L], 3L)
  expect_match(headers, "#define CBW_DEGREE_SEARCHI 33", fixed = TRUE)
  expect_match(headers, "#define CBW_CATCOMPI 34", fixed = TRUE)
  expect_match(headers, "#define CDBW_CATCOMPI 32", fixed = TRUE)
  expect_match(
    ingress,
    paste0(
      "np_conditional_density_bw_categorical_compress_extern =\n",
      "      myopti[CBW_CATCOMPI];"
    ),
    fixed = TRUE
  )
  expect_match(
    ingress,
    paste0(
      "np_conditional_distribution_bw_categorical_compress_extern =\n",
      "      myopti[CDBW_CATCOMPI];"
    ),
    fixed = TRUE
  )
  expect_match(
    ingress,
    "C_np_density_conditional_bw: categorical compression must be TRUE or FALSE",
    fixed = TRUE
  )
  expect_match(
    ingress,
    "C_np_distribution_conditional_bw: categorical compression must be TRUE or FALSE",
    fixed = TRUE
  )
})

test_that("conditional beta bandwidth ingress is open for both estimators", {
  root <- locate_beta_conditional_bandwidth_sources()
  skip_if(is.null(root), "package sources unavailable")
  ingress <- paste(
    readLines(file.path(root, "src", "np.c"), warn = FALSE),
    collapse = "\n"
  )

  expect_false(grepl(
    "C_np_density_conditional_bw: beta bandwidth selection supports only local-constant fitting",
    ingress,
    fixed = TRUE
  ))
  expect_false(grepl(
    "C_np_density_conditional_bw: beta bandwidth selection requires continuous X and Y variables only",
    ingress,
    fixed = TRUE
  ))
  expect_false(grepl(
    "C_np_distribution_conditional_bw: beta bandwidth selection supports only local-constant fitting",
    ingress,
    fixed = TRUE
  ))
  expect_false(grepl(
    "C_np_distribution_conditional_bw: beta bandwidth selection requires continuous X and Y variables only",
    ingress,
    fixed = TRUE
  ))

  density_starts <- gregexpr(
    "void np_density_conditional_bw(double * c_uno",
    ingress,
    fixed = TRUE
  )[[1L]]
  distribution_starts <- gregexpr(
    "void np_distribution_conditional_bw(double * c_uno",
    ingress,
    fixed = TRUE
  )[[1L]]
  density_start <- tail(density_starts[density_starts > 0L], 1L)
  distribution_start <- tail(
    distribution_starts[distribution_starts > 0L],
    1L
  )
  expect_gt(density_start, 0L)
  expect_gt(distribution_start, density_start)
  density_ingress <- substr(
    ingress,
    density_start,
    distribution_start - 1L
  )
  allow_hits <- gregexpr(
    "NP_BETA_BW_ALLOW_CATEGORICAL",
    density_ingress,
    fixed = TRUE
  )[[1L]]
  expect_length(allow_hits[allow_hits > 0L], 2L)
  expect_false(grepl(
    "NP_BETA_BW_CONTINUOUS_ONLY",
    density_ingress,
    fixed = TRUE
  ))
  distribution_ingress <- substr(
    ingress,
    distribution_start,
    nchar(ingress)
  )
  distribution_allow_hits <- gregexpr(
    "NP_BETA_BW_ALLOW_CATEGORICAL",
    distribution_ingress,
    fixed = TRUE
  )[[1L]]
  expect_gte(length(distribution_allow_hits[
    distribution_allow_hits > 0L
  ]), 2L)
  expect_match(
    paste(readLines(file.path(root, "R", "condbandwidth.R"),
                    warn = FALSE), collapse = "\n"),
    "allow.automatic.full = TRUE",
    fixed = TRUE
  )
})

test_that("distribution degree search owns automatic beta placeholders", {
  root <- locate_beta_conditional_bandwidth_sources()
  skip_if(is.null(root), "package sources unavailable")
  distribution_r <- paste(
    readLines(file.path(root, "R", "np.condistribution.bw.R"),
              warn = FALSE),
    collapse = "\n"
  )
  search_start <- regexpr(
    ".npcdistbw_nomad_search <- function",
    distribution_r,
    fixed = TRUE
  )[[1L]]
  expect_gt(search_start, 0L)
  search_text <- substring(distribution_r, search_start)
  expect_match(
    search_text,
    paste0(
      "bandwidth.compute = identical(reg.args$cxkertype, \"beta\") ||\n",
      "      identical(reg.args$cykertype, \"beta\")"
    ),
    fixed = TRUE
  )
})

test_that("conditional bandwidth objectives are neutral to dormant compression state", {
  .ensure_beta_conditional_ingress_pool()
  old <- options(np.messages = FALSE, np.categorical.compress = FALSE)
  on.exit(options(old), add = TRUE)

  x <- data.frame(x = c(0.03, 0.09, 0.16, 0.27, 0.39, 0.52,
                        0.64, 0.73, 0.82, 0.9, 0.96))
  y <- data.frame(y = c(0.04, 0.12, 0.21, 0.33, 0.47, 0.58,
                        0.69, 0.78, 0.86, 0.93, 0.98))
  density_args <- list(
    xdat = x, ydat = y, bws = c(0.15, 0.17),
    bandwidth.compute = FALSE, bwmethod = "cv.ml",
    cxkertype = "beta", cxkerorder = 4L,
    cxkerbound = "fixed", cxkerlb = 0, cxkerub = 1,
    cykertype = "beta", cykerorder = 6L,
    cykerbound = "fixed", cykerlb = 0, cykerub = 1
  )
  distribution_args <- density_args
  distribution_args$bwmethod <- NULL

  options(np.categorical.compress = FALSE)
  density_bw <- do.call(npcdensbw, density_args)
  distribution_bw <- do.call(npcdistbw, distribution_args)
  density_dense <- npRmpi:::.npcdensbw_eval_only(x, y, density_bw)
  distribution_dense <- npRmpi:::.npcdistbw_eval_only(x, y, bws = distribution_bw)
  options(np.categorical.compress = TRUE)
  density_compressed <- npRmpi:::.npcdensbw_eval_only(x, y, density_bw)
  distribution_compressed <- npRmpi:::.npcdistbw_eval_only(
    x, y, bws = distribution_bw
  )

  expect_identical(density_compressed$objective, density_dense$objective)
  expect_identical(
    distribution_compressed$objective, distribution_dense$objective
  )
})

test_that("conditional bandwidth ingress rejects malformed compression state", {
  .ensure_beta_conditional_ingress_pool()
  old <- options(np.messages = FALSE, np.categorical.compress = TRUE)
  on.exit(options(old), add = TRUE)

  x <- data.frame(x = seq(0.05, 0.95, length.out = 9L))
  y <- data.frame(y = seq(0.08, 0.92, length.out = 9L))
  args <- list(
    xdat = x, ydat = y, bws = c(0.16, 0.18),
    bandwidth.compute = FALSE,
    cxkertype = "beta", cxkerbound = "fixed",
    cxkerlb = 0, cxkerub = 1,
    cykertype = "beta", cykerbound = "fixed",
    cykerlb = 0, cykerub = 1
  )
  message <- paste0(
    "option 'np.categorical.compress' must be a single ",
    "non-missing logical value"
  )

  density_bw <- do.call(npcdensbw, args)
  distribution_bw <- do.call(npcdistbw, args)
  options(np.categorical.compress = NA)

  expect_error(
    npRmpi:::.npcdensbw_eval_only(x, y, density_bw), message, fixed = TRUE
  )
  expect_error(
    npRmpi:::.npcdistbw_eval_only(x, y, bws = distribution_bw),
    message,
    fixed = TRUE
  )
})
