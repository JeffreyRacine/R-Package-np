locate_demo_parser <- function() {
  candidates <- c(
    test_path("..", "..", "inst", "demo_tools", "parse_demo_results.R"),
    test_path("..", "..", "..", "inst", "demo_tools", "parse_demo_results.R"),
    system.file("demo_tools", "parse_demo_results.R", package = "npRmpi")
  )
  hits <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (!length(hits)) return("")
  normalizePath(hits[[1L]], mustWork = TRUE)
}

demo_compute_specs <- function() {
  list(
    list(path = c("serial"), mode = "serial", extra = ""),
    list(path = c("session", "slaves_01"), mode = "session", extra = "slaves=1"),
    list(path = c("session", "slaves_02"), mode = "session", extra = "slaves=2"),
    list(path = c("session", "slaves_03"), mode = "session", extra = "slaves=3"),
    list(path = c("mpi_launch", "ranks_02", "attach"), mode = "attach", extra = "ranks=2"),
    list(path = c("mpi_launch", "ranks_03", "attach"), mode = "attach", extra = "ranks=3"),
    list(path = c("mpi_launch", "ranks_04", "attach"), mode = "attach", extra = "ranks=4"),
    list(path = c("mpi_launch", "ranks_02", "profile"), mode = "profile", extra = "ranks=2"),
    list(path = c("mpi_launch", "ranks_03", "profile"), mode = "profile", extra = "ranks=3"),
    list(path = c("mpi_launch", "ranks_04", "profile"), mode = "profile", extra = "ranks=4")
  )
}

write_demo_fixture <- function(root, case, descriptors = "", offset = 0) {
  specs <- demo_compute_specs()
  for (i in seq_along(specs)) {
    spec <- specs[[i]]
    directory <- do.call(file.path, c(list(root), as.list(spec$path)))
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    line <- paste(
      "DEMO_RESULT",
      paste0("demo=", case),
      paste0("mode=", spec$mode),
      spec$extra,
      "n=100",
      "default_n=100",
      sprintf("elapsed=%.3f", offset + i / 10),
      "family=fixture",
      paste0("case=", case),
      "tier=sentinel",
      descriptors
    )
    writeLines(line, file.path(directory, paste0(case, ".Rout")))
  }
}

run_demo_parser <- function(root, output) {
  parser <- locate_demo_parser()
  if (!nzchar(parser)) return(list(status = NA_integer_, output = ""))
  result <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c("--vanilla", parser, root, output),
    stdout = TRUE,
    stderr = TRUE
  ))
  list(status = if (is.null(attr(result, "status"))) 0L else attr(result, "status"),
       output = result)
}

test_that("wide demo timing retains cases with missing optional descriptors", {
  skip_if(!nzchar(locate_demo_parser()), "demo timing parser unavailable")
  root <- tempfile("demo-wide-valid-")
  output <- tempfile("demo-wide-output-")
  dir.create(root)
  dir.create(output)

  write_demo_fixture(
    root,
    "complete_case",
    "regtype=ll bwmethod=cv.ls nomad=FALSE degree=NA degree.max=NA selected.degree=1 bwtype=fixed",
    offset = 0
  )
  write_demo_fixture(root, "missing_case", offset = 10)

  parsed <- run_demo_parser(root, output)
  expect_identical(parsed$status, 0L, info = paste(parsed$output, collapse = "\n"))

  long <- read.csv(file.path(output, "demo_results.csv"), stringsAsFactors = FALSE)
  wide <- read.csv(file.path(output, "demo_results_wide.csv"), stringsAsFactors = FALSE)
  expect_equal(nrow(long), 20L)
  expect_equal(nrow(wide), 2L)
  expect_setequal(wide$case, c("complete_case", "missing_case"))

  timing_columns <- c(
    "seconds_serial",
    sprintf("seconds_session_s%02d", 1:3),
    sprintf("seconds_attach_r%02d", 2:4),
    sprintf("seconds_profile_r%02d", 2:4)
  )
  expect_true(all(timing_columns %in% names(wide)))

  for (i in seq_len(nrow(long))) {
    compute <- if (long$mode[[i]] == "serial") {
      "serial"
    } else if (long$mode[[i]] == "session") {
      sprintf("session_s%02d", long$slaves[[i]])
    } else {
      sprintf("%s_r%02d", long$mode[[i]], long$ranks[[i]])
    }
    row <- match(long$case[[i]], wide$case)
    expect_equal(wide[[paste0("seconds_", compute)]][[row]], long$elapsed[[i]],
                 tolerance = 0)
  }

  missing_row <- match("missing_case", wide$case)
  expect_true(is.na(wide$regtype[[missing_row]]))
  expect_true(is.na(wide$bwmethod[[missing_row]]))
  expect_true(is.na(wide$selected.degree[[missing_row]]))
})

test_that("wide demo timing rejects conflicting descriptors", {
  skip_if(!nzchar(locate_demo_parser()), "demo timing parser unavailable")
  root <- tempfile("demo-wide-conflict-")
  output <- tempfile("demo-wide-conflict-output-")
  dir.create(root)
  dir.create(output)

  write_demo_fixture(root, "conflict_case", "regtype=ll bwmethod=cv.ls", offset = 0)
  serial <- file.path(root, "serial", "conflict_case.Rout")
  line <- readLines(serial, warn = FALSE)
  writeLines(sub("regtype=ll", "regtype=lc", line, fixed = TRUE), serial)

  parsed <- run_demo_parser(root, output)
  expect_false(is.na(parsed$status))
  expect_false(identical(parsed$status, 0L))
  expect_match(paste(parsed$output, collapse = "\n"),
               "conflicting wide descriptor regtype for case=conflict_case",
               fixed = TRUE)
})

test_that("wide demo timing rejects duplicate compute cells", {
  skip_if(!nzchar(locate_demo_parser()), "demo timing parser unavailable")
  root <- tempfile("demo-wide-duplicate-")
  output <- tempfile("demo-wide-duplicate-output-")
  dir.create(root)
  dir.create(output)

  write_demo_fixture(root, "duplicate_case", offset = 0)
  file.copy(file.path(root, "serial", "duplicate_case.Rout"),
            file.path(root, "serial", "duplicate_case_copy.Rout"))

  parsed <- run_demo_parser(root, output)
  expect_false(is.na(parsed$status))
  expect_false(identical(parsed$status, 0L))
  expect_match(paste(parsed$output, collapse = "\n"),
               "duplicate wide timing cell for case=duplicate_case compute=serial",
               fixed = TRUE)
})
