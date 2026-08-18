locate_exact_field_source <- function(file) {
  roots <- unique(c(
    Sys.getenv("R_PACKAGE_SOURCE", ""),
    testthat::test_path("..", ".."),
    testthat::test_path("..", "..", ".."),
    getwd(),
    file.path(getwd(), "..")
  ))
  paths <- file.path(roots[nzchar(roots)], "R", file)
  paths <- paths[file.exists(paths)]
  if (!length(paths))
    return(NULL)
  paths[[1L]]
}

test_that("campaign-added public object reads use exact extraction", {
  expected <- list(
    "np.condensity.bw.R" = c(
      'bws[["cykerorder", exact = TRUE]]',
      'bws[["xbw", exact = TRUE]][bws[["ixcon", exact = TRUE]]]',
      'bws[["ybw", exact = TRUE]][bws[["iycon", exact = TRUE]]]',
      'bws[["ybw", exact = TRUE]][bws[["iyuno", exact = TRUE]]]',
      'bws[["ybw", exact = TRUE]][bws[["iyord", exact = TRUE]]]',
      'bws[["xbw", exact = TRUE]][bws[["ixuno", exact = TRUE]]]',
      'bws[["xbw", exact = TRUE]][bws[["ixord", exact = TRUE]]]'
    ),
    "np.condistribution.bw.R" = c(
      'bws[["cykerorder", exact = TRUE]]',
      'bws[["xbw", exact = TRUE]][bws[["ixcon", exact = TRUE]]]',
      'bws[["ybw", exact = TRUE]][bws[["iycon", exact = TRUE]]]',
      'bws[["ybw", exact = TRUE]][bws[["iyuno", exact = TRUE]]]',
      'bws[["ybw", exact = TRUE]][bws[["iyord", exact = TRUE]]]',
      'bws[["xbw", exact = TRUE]][bws[["ixuno", exact = TRUE]]]',
      'bws[["xbw", exact = TRUE]][bws[["ixord", exact = TRUE]]]'
    ),
    "np.kernel.R" = c(
      'bws[["nuno", exact = TRUE]] + bws[["nord", exact = TRUE]]'
    ),
    "np.plot.helpers.R" = c(
      'fit[[if (cdf) "condist" else "condens", exact = TRUE]]',
      'bws[["type", exact = TRUE]]'
    ),
    "np.reghat.R" = c(
      'bws[["icon", exact = TRUE]]',
      'bws[["ckertype", exact = TRUE]]',
      'bws[["ckerorder", exact = TRUE]]',
      'bws[["ncon", exact = TRUE]]',
      'bws[["nuno", exact = TRUE]] + bws[["nord", exact = TRUE]]',
      'bws[["ckerlb", exact = TRUE]][bws[["icon", exact = TRUE]]]',
      'bws[["ckerub", exact = TRUE]][bws[["icon", exact = TRUE]]]'
    ),
    "np.regression.R" = c(
      'bws[["nuno", exact = TRUE]]',
      'bws[["nord", exact = TRUE]]'
    )
  )

  for (file in names(expected)) {
    path <- locate_exact_field_source(file)
    skip_if(is.null(path), paste("package source unavailable for", file))
    source <- paste(readLines(path, warn = FALSE), collapse = "\n")
    for (literal in expected[[file]])
      expect_match(source, literal, fixed = TRUE, info = file)
  }
})

test_that("a similarly prefixed bandwidth field cannot activate reghat", {
  package <- if ("npRmpi" %in% loadedNamespaces()) "npRmpi" else "np"
  candidate <- getFromNamespace(".npreghat_native_apply_candidate", package)
  arguments <- list(
    output = "apply",
    y = matrix(1, nrow = 3L, ncol = 2L),
    regtype.engine = "lp",
    degree = 1L,
    basis = "glp",
    bernstein.basis = FALSE,
    s = 0L,
    leave.one.out = FALSE
  )

  prefix.only <- c(list(bws = list(icon.shadow = TRUE, ckertype = "beta")),
                   arguments)
  exact <- c(list(bws = list(icon = TRUE, ckertype = "beta")), arguments)

  expect_false(do.call(candidate, prefix.only))
  expect_true(do.call(candidate, exact))
})
