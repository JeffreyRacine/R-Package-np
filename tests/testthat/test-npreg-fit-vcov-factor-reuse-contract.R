locate_lp_source <- function(name) {
  candidates <- c(
    test_path("..", "..", "src", name),
    test_path("..", "..", "..", "src", name),
    file.path(Sys.getenv("R_PACKAGE_DIR", ""), "src", name),
    file.path(Sys.getenv("R_PACKAGE_SOURCE", ""), "src", name),
    file.path(getwd(), "src", name),
    file.path(getwd(), "..", "src", name)
  )
  candidates <- unique(candidates[nzchar(candidates)])
  hits <- candidates[file.exists(candidates)]
  if (length(hits) == 0L) NULL else hits[[1L]]
}

function_body <- function(lines, signature) {
  start <- grep(signature, lines)
  expect_length(start, 1L)
  stops <- grep("^}$", lines)
  stop <- stops[stops > start][1L]
  expect_length(stop, 1L)
  paste(lines[start:stop], collapse = "\n")
}

test_that("LP fit covariance reuses only a validated retained factorization", {
  jksum_file <- locate_lp_source("jksum.c")
  solve_file <- locate_lp_source("jksum_lp_solve.c")
  header_file <- locate_lp_source("jksum_lp_solve.h")
  skip_if(
    any(vapply(
      list(jksum_file, solve_file, header_file),
      is.null,
      logical(1L)
    )),
    "LP source files unavailable in this test context"
  )

  jksum_lines <- readLines(jksum_file, warn = FALSE)
  solve_lines <- readLines(solve_file, warn = FALSE)
  header_lines <- readLines(header_file, warn = FALSE)

  expect_equal(
    sum(grepl("np_lp_solve_workspace_solve_factored\\(", jksum_lines)),
    1L
  )
  expect_true(any(grepl("int factor_ready;", header_lines, fixed = TRUE)))
  expect_true(any(grepl("int factor_p;", header_lines, fixed = TRUE)))

  solve_body <- function_body(
    solve_lines,
    "^int np_lp_solve_workspace_solve\\("
  )
  expect_true(grepl("workspace->factor_ready = 0;", solve_body, fixed = TRUE))
  expect_true(grepl("workspace->factor_ready = 1;", solve_body, fixed = TRUE))
  expect_true(grepl("workspace->factor_p = p;", solve_body, fixed = TRUE))

  factored_body <- function_body(
    solve_lines,
    "^int np_lp_solve_workspace_solve_factored\\("
  )
  expect_true(grepl("!workspace->factor_ready", factored_body, fixed = TRUE))
  expect_true(grepl("workspace->factor_p != p", factored_body, fixed = TRUE))
  expect_true(grepl("F77_CALL(dgetrs)", factored_body, fixed = TRUE))
  expect_false(grepl("F77_CALL(dgesv)", factored_body, fixed = TRUE))
  expect_false(grepl("F77_CALL(dgetrf)", factored_body, fixed = TRUE))
  expect_false(grepl("gram_source", factored_body, fixed = TRUE))
})
