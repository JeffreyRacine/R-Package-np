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
    header = file.path(roots[[1L]], "src", "headers.h")
  )
}

test_that("LP delete-row denominators use the shared signed guard", {
  paths <- locate_signed_delete_sources()
  skip_if(is.null(paths), "package C sources unavailable in this test context")

  header <- paste(readLines(paths$header, warn = FALSE), collapse = "\n")

  expect_match(header, "static inline double NZD(const double a)", fixed = TRUE)
  expect_match(
    header,
    "return (a < 0.0) ? MIN(-DBL_EPSILON, a) : MAX(DBL_EPSILON, a);",
    fixed = TRUE
  )
})

test_that("promoted conditional CVLS paired rows use signed deletion", {
  paths <- locate_signed_delete_sources()
  skip_if(is.null(paths), "package C sources unavailable in this test context")

  source <- paste(readLines(paths$jksum, warn = FALSE), collapse = "\n")
  expect_equal(lengths(regmatches(
    source,
    gregexpr(
      "NZD\\(1\\.0 - full_rows_out\\[i\\]\\[eval_idx\\]\\)",
      source,
      perl = TRUE
    )
  )), 2L)
  expect_false(grepl(
    "NZD_POS(1.0 - full_rows_out[i][eval_idx])",
    source,
    fixed = TRUE
  ))
})

test_that("higher-order fixed and generalized-NN CVLS match signed-WLS oracles", {
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
      oracle = 0.8084230657971182
    ),
    list(
      bwtype = "fixed", kernel = "gaussian", order = 4L,
      bernstein = FALSE, degree = c(2L, 2L),
      bws = c(0.28, 0.31, 0.27),
      oracle = 1.1317864794839687
    ),
    list(
      bwtype = "generalized_nn", kernel = "gaussian", order = 8L,
      bernstein = TRUE, degree = c(2L, 2L),
      bws = c(31L, 34L, 29L),
      oracle = 2.4431117127599182
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
