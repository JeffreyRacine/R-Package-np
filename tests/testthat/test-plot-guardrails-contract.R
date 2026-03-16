test_that("plot helper runtime guard requires non-null bws", {
  guard <- getFromNamespace(".np_plot_require_bws", "npRmpi")
  expect_error(
    guard(NULL, "unit-test"),
    "required argument 'bws' is missing or NULL"
  )
  expect_invisible(guard(list(type = "fixed"), "unit-test"))
})

test_that("plot runtime files avoid forbidden *bw( calls", {
  root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  files <- c(
    file.path(root, "R", "np.plot.helpers.R"),
    Sys.glob(file.path(root, "R", "np.plot.engine*.R"))
  )
  files <- unique(files[file.exists(files)])
  skip_if(length(files) == 0L, "source R files unavailable in installed test context")

  forbidden <- "\\b[A-Za-z0-9._]+bw\\s*\\("
  allowed <- c(
    "\\.np_indexhat_rbw\\s*\\(",
    "\\.np_indexhat_kbw\\s*\\(",
    "\\.npcdhat_make_xbw\\s*\\(",
    "\\.npcdhat_make_xkbw\\s*\\(",
    "\\.npcdhat_make_ybw\\s*\\("
  )
  offenders <- character()

  for (f in files) {
    raw <- readLines(f, warn = FALSE)
    code <- sub("#.*$", "", raw)
    hit_idx <- which(grepl(forbidden, code, perl = TRUE))
    if (!length(hit_idx)) next

    for (i in hit_idx) {
      ok <- any(vapply(allowed, function(p) grepl(p, code[[i]], perl = TRUE), logical(1)))
      if (!ok) {
        offenders <- c(offenders, sprintf("%s:%d: %s", basename(f), i, trimws(raw[[i]])))
      }
    }
  }

  expect_equal(length(offenders), 0, info = paste(offenders, collapse = "\n"))
})

test_that("npRmpi plot runtime files avoid serial namespace-bridge calls", {
  root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  files <- c(
    file.path(root, "R", "np.plot.helpers.R"),
    Sys.glob(file.path(root, "R", "np.plot*.R"))
  )
  files <- unique(files[file.exists(files)])
  skip_if(length(files) == 0L, "source R files unavailable in installed test context")

  offenders <- character()
  for (f in files) {
    raw <- readLines(f, warn = FALSE)
    bridge.pat <- paste0("np", "::")
    idx <- which(grepl(bridge.pat, sub("#.*$", "", raw), fixed = TRUE))
    if (length(idx)) {
      offenders <- c(offenders, sprintf("%s:%d: %s", basename(f), idx, trimws(raw[idx])))
    }
  }

  expect_equal(length(offenders), 0, info = paste(offenders, collapse = "\n"))
})

test_that("npRmpi plot runtime files avoid silent remap/downgrade patterns", {
  root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  files <- c(
    file.path(root, "R", "np.plot.helpers.R"),
    Sys.glob(file.path(root, "R", "np.plot.engine*.R"))
  )
  files <- unique(files[file.exists(files)])
  skip_if(length(files) == 0L, "source R files unavailable in installed test context")

  offenders <- character()
  pat.assign <- "plot\\.errors\\.method\\s*=\\s*\"none\""
  pat.warn <- "Proceeding without"

  for (f in files) {
    raw <- readLines(f, warn = FALSE)
    code <- sub("#.*$", "", raw)
    idx <- which(grepl(pat.assign, code, perl = TRUE) | grepl(pat.warn, raw, fixed = TRUE))
    if (length(idx)) {
      offenders <- c(offenders, sprintf("%s:%d: %s", basename(f), idx, trimws(raw[idx])))
    }
  }

  expect_equal(length(offenders), 0, info = paste(offenders, collapse = "\n"))
})

test_that("npRmpi runtime and demos avoid nested mpi.bcast.cmd(plot(...)) calls", {
  root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  files <- c(
    Sys.glob(file.path(root, "R", "*.R")),
    Sys.glob(file.path(root, "demo", "*.R"))
  )
  files <- unique(files[file.exists(files)])
  skip_if(length(files) == 0L, "source R/demo files unavailable in installed test context")

  pat <- "mpi\\.bcast\\.cmd\\s*\\(\\s*(?:[A-Za-z0-9._]+::)?plot\\s*\\("
  offenders <- character()

  for (f in files) {
    raw <- readLines(f, warn = FALSE)
    code <- paste(sub("#.*$", "", raw), collapse = "\n")
    if (grepl(pat, code, perl = TRUE))
      offenders <- c(offenders, basename(f))
  }

  expect_equal(length(offenders), 0, info = paste(offenders, collapse = ", "))
})
