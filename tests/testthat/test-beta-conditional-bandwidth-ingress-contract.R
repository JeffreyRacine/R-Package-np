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

test_that("conditional bandwidth objectives are neutral to dormant compression state", {
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
  density_dense <- np:::.npcdensbw_eval_only(x, y, density_bw)
  distribution_dense <- np:::.npcdistbw_eval_only(x, y, bws = distribution_bw)
  options(np.categorical.compress = TRUE)
  density_compressed <- np:::.npcdensbw_eval_only(x, y, density_bw)
  distribution_compressed <- np:::.npcdistbw_eval_only(
    x, y, bws = distribution_bw
  )

  expect_identical(density_compressed$objective, density_dense$objective)
  expect_identical(
    distribution_compressed$objective, distribution_dense$objective
  )
})

test_that("conditional bandwidth ingress rejects malformed compression state", {
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
    np:::.npcdensbw_eval_only(x, y, density_bw), message, fixed = TRUE
  )
  expect_error(
    np:::.npcdistbw_eval_only(x, y, bws = distribution_bw),
    message,
    fixed = TRUE
  )
})
