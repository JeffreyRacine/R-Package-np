locate_signed_delete_sources <- function() {
  roots <- c(
    test_path("..", ".."),
    test_path("..", "..", ".."),
    Sys.getenv("R_PACKAGE_DIR", ""),
    Sys.getenv("R_PACKAGE_SOURCE", ""),
    getwd(),
    file.path(getwd(), "..")
  )
  roots <- unique(roots[nzchar(roots)])
  roots <- roots[file.exists(file.path(roots, "src", "jksum.c"))]
  if (!length(roots)) return(NULL)
  list(
    jksum = file.path(roots[[1L]], "src", "jksum.c"),
    row = file.path(roots[[1L]], "src", "jksum_lp_row.c"),
    row_header = file.path(roots[[1L]], "src", "jksum_lp_row.h")
  )
}

test_that("LP delete-row denominators retain every finite nonzero signed value", {
  paths <- locate_signed_delete_sources()
  skip_if(is.null(paths), "package C sources unavailable in this test context")

  header <- paste(readLines(paths$row_header, warn = FALSE), collapse = "\n")
  row <- paste(readLines(paths$row, warn = FALSE), collapse = "\n")

  expect_match(header, "static inline int np_lp_delete_denominator", fixed = TRUE)
  expect_match(header, "*denominator = 1.0 - leverage;", fixed = TRUE)
  expect_match(header, "isfinite(*denominator)", fixed = TRUE)
  expect_match(header, "(*denominator != 0.0)", fixed = TRUE)
  expect_false(grepl("DBL_EPSILON", header, fixed = TRUE))
  expect_false(grepl("np_lp_delete_smoother_row", row, fixed = TRUE))
  expect_false(grepl("np_lp_delete_smoother_row", header, fixed = TRUE))
})

test_that("canonical conditional CVLS rows use exact signed deletion", {
  paths <- locate_signed_delete_sources()
  skip_if(is.null(paths), "package C sources unavailable in this test context")

  source <- paste(readLines(paths$jksum, warn = FALSE), collapse = "\n")
  expect_match(
    source,
    "np_lp_delete_denominator(rows_out[i][eval_idx], &den)",
    fixed = TRUE
  )
  expect_false(grepl("np_conditional_x_weight_block_pair", source,
                     fixed = TRUE))
  expect_false(grepl("np_lp_delete_smoother_row", source, fixed = TRUE))
  expect_false(grepl(
    "NZD(1.0 - full_rows_out[i][eval_idx])",
    source,
    fixed = TRUE
  ))
  expect_false(grepl(
    "NZD_POS(1.0 - full_rows_out[i][eval_idx])",
    source,
    fixed = TRUE
  ))
})

test_that("higher-order fixed and generalized-NN CVLS retain qualified transcripts", {
  skip_if_not_installed("np")
  suppressPackageStartupMessages(library(np))
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  set.seed(2026072701L)
  n <- 96L
  xdat <- data.frame(
    x1 = runif(n, 0.05, 0.95),
    x2 = runif(n, 0.05, 0.95)
  )
  ydat <- data.frame(
    y = 0.5 + 0.22 * sin(2 * pi * xdat$x1 - 0.7 * xdat$x2) +
      rnorm(n, sd = 0.08)
  )
  cases <- list(
    list(
      bwtype = "fixed", kernel = "epanechnikov", order = 6L,
      bernstein = FALSE, degree = c(1L, 1L),
      bws = c(0.52, 0.56, 0.48),
      oracle = 0.78513338511644237
    ),
    list(
      bwtype = "fixed", kernel = "gaussian", order = 4L,
      bernstein = FALSE, degree = c(2L, 2L),
      bws = c(0.28, 0.31, 0.27),
      oracle = -0.13149110780478809
    ),
    list(
      bwtype = "generalized_nn", kernel = "gaussian", order = 8L,
      bernstein = TRUE, degree = c(2L, 2L),
      bws = c(31L, 34L, 29L),
      oracle = 0.79852466654749321
    )
  )

  for (case in cases) {
    bws <- npcdensbw(
      xdat = xdat,
      ydat = ydat,
      bws = case$bws,
      bandwidth.compute = FALSE,
      bwmethod = "cv.ls",
      bwtype = case$bwtype,
      regtype = "lp",
      degree = case$degree,
      bernstein.basis = case$bernstein,
      cxkertype = case$kernel,
      cxkerorder = case$order,
      cykertype = "gaussian",
      cykerorder = 2L
    )
    objective <- getFromNamespace(".npcdensbw_eval_only", "np")(
      xdat = xdat,
      ydat = ydat,
      bws = bws
    )$objective
    expect_true(is.finite(objective))
    expect_equal(as.numeric(objective), case$oracle, tolerance = 2e-8)
  }
})

test_that("remaining conditional higher-order LOO routes retain qualified transcripts", {
  skip_if_not_installed("np")
  suppressPackageStartupMessages(library(np))
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  set.seed(2026072702L)
  n <- 72L
  xdat <- data.frame(
    x1 = runif(n, 0.05, 0.95),
    x2 = runif(n, 0.05, 0.95)
  )
  ydat <- data.frame(
    y = 0.4 + 0.3 * sin(2 * pi * xdat$x1) * cos(pi * xdat$x2) +
      rnorm(n, sd = 0.09)
  )
  # These are adjacency transcripts for the high-order signed-row owners.
  # Exact donor-fold NN geometry is proved independently by the dedicated
  # adaptive/generalized conditional CV oracle contracts.
  cases <- list(
    list(
      route = "cdens", method = "cv.ml", bwtype = "fixed",
      kernel = "gaussian", order = 4L, bernstein = FALSE,
      degree = c(1L, 1L), bws = c(0.25, 0.29, 0.27),
      oracle = -1391.733240027958
    ),
    list(
      route = "cdens", method = "cv.ls", bwtype = "adaptive_nn",
      kernel = "epanechnikov", order = 8L, bernstein = TRUE,
      degree = c(2L, 2L), bws = c(25L, 29L, 27L),
      oracle = -46.797292511161814
    ),
    list(
      route = "cdist", method = "cv.ls", bwtype = "generalized_nn",
      kernel = "gaussian", order = 8L, bernstein = FALSE,
      degree = c(1L, 1L), bws = c(25L, 29L, 27L),
      oracle = 3.066587903283053
    )
  )

  for (case in cases) {
    constructor <- if (identical(case$route, "cdens")) npcdensbw else npcdistbw
    bw <- constructor(
      xdat = xdat,
      ydat = ydat,
      bws = case$bws,
      bandwidth.compute = FALSE,
      bwmethod = case$method,
      bwtype = case$bwtype,
      regtype = "lp",
      basis = "glp",
      degree = case$degree,
      bernstein.basis = case$bernstein,
      cxkertype = case$kernel,
      cxkerorder = case$order,
      cykertype = "gaussian",
      cykerorder = 2L
    )
    objective <- if (identical(case$route, "cdens")) {
      np:::.npcdensbw_eval_only(xdat, ydat, bw)$objective
    } else {
      np:::.npcdistbw_eval_only(
        xdat, ydat, bws = bw, do.full.integral = TRUE
      )$objective
    }
    expect_true(is.finite(objective))
    expect_equal(as.numeric(objective), case$oracle, tolerance = 5e-8)
  }
})
