locate_mpi_lp_source <- function() {
  candidates <- c(
    test_path("..", "..", "src", "jksum.c"),
    test_path("..", "..", "..", "src", "jksum.c"),
    file.path(Sys.getenv("R_PACKAGE_DIR", ""), "src", "jksum.c"),
    file.path(Sys.getenv("R_PACKAGE_SOURCE", ""), "src", "jksum.c"),
    file.path(getwd(), "src", "jksum.c"),
    file.path(getwd(), "..", "src", "jksum.c")
  )
  candidates <- unique(candidates[nzchar(candidates)])
  hits <- candidates[file.exists(candidates)]
  if (!length(hits))
    return(NULL)
  hits[[1L]]
}

test_that("MPI directed width-one rows use implicit scalar accumulation", {
  source_file <- locate_mpi_lp_source()
  skip_if(is.null(source_file), "source file src/jksum.c unavailable")
  source <- paste(readLines(source_file, warn = FALSE), collapse = "\n")

  owner_start <- regexpr(
    paste0(
      "if(use_mpi_transport){\n",
      "      double * const sj = moments_acc"
    ),
    source,
    fixed = TRUE
  )
  expect_gt(owner_start, 0L)
  owner_source <- substr(source, owner_start, nchar(source))
  scalar_start <- regexpr("if(nterms == 1){", owner_source, fixed = TRUE)
  wider_start <- regexpr(
    paste0(
      "} else {\n",
      "        const NPLPOwnedRowContext owned_context"
    ),
    owner_source,
    fixed = TRUE
  )
  expect_gt(scalar_start, 0L)
  expect_gt(wider_start, scalar_start)
  scalar <- substr(owner_source, scalar_start, wider_start - 1L)

  expect_true(grepl("double sj0 = sj[0]", scalar, fixed = TRUE))
  expect_true(grepl("double tj0 = tj[0]", scalar, fixed = TRUE))
  expect_true(grepl("tj0 += w*vector_Y[ii]", scalar, fixed = TRUE))
  expect_true(grepl("sj0 += w", scalar, fixed = TRUE))
  expect_false(grepl("(^|[^[:alnum:]_])basis\\[", scalar, perl = TRUE))
  expect_false(grepl("F77_", scalar, fixed = TRUE))
})

test_that("MPI sparse-tree directed pairs short-circuit width one", {
  source_file <- locate_mpi_lp_source()
  skip_if(is.null(source_file), "source file src/jksum.c unavailable")
  source <- paste(readLines(source_file, warn = FALSE), collapse = "\n")

  start <- regexpr("static inline void np_lp_accumulate_row(", source,
                   fixed = TRUE)
  stop <- regexpr("static int np_lp_tree_oracle_enabled(", source,
                  fixed = TRUE)
  expect_gt(start, 0L)
  expect_gt(stop, start)
  owner <- substr(source, start, stop - 1L)
  scalar_start <- regexpr("if(nterms == 1)", owner, fixed = TRUE)
  wider_start <- regexpr("for(a = 0; a < nterms; a++)", owner, fixed = TRUE)
  expect_gt(scalar_start, 0L)
  expect_gt(wider_start, scalar_start)
  scalar <- substr(owner, scalar_start, wider_start - 1L)

  expect_true(grepl("tj[0] += w*yi", scalar, fixed = TRUE))
  expect_true(grepl("sj[0] += w", scalar, fixed = TRUE))
  expect_false(grepl("(^|[^[:alnum:]_])basis\\[", scalar, perl = TRUE))
  expect_false(grepl("F77_", scalar, fixed = TRUE))
})
