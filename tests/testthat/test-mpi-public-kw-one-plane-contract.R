locate_mpi_public_kw_source <- function() {
  roots <- unique(c(
    testthat::test_path("..", ".."),
    testthat::test_path("..", "..", ".."),
    Sys.getenv("R_PACKAGE_DIR", ""),
    Sys.getenv("R_PACKAGE_SOURCE", ""),
    getwd(),
    file.path(getwd(), "..")
  ))
  paths <- file.path(roots[nzchar(roots)], "src", "jksum.c")
  paths <- paths[file.exists(paths)]
  if (!length(paths))
    return(NULL)
  paths[[1L]]
}

mpi_public_kw_source_region <- function(lines, start_pattern, stop_pattern) {
  start <- grep(start_pattern, lines)
  if (length(start) != 1L)
    stop("expected exactly one kernel-weight owner start marker")
  stop <- grep(stop_pattern, lines)
  stop <- stop[stop > start][1L]
  if (length(stop) != 1L || is.na(stop))
    stop("kernel-weight owner stop marker is missing")
  paste(lines[start:(stop - 1L)], collapse = "\n")
}

test_that("MPI public kernel weights use one plane and typed collectives", {
  path <- locate_mpi_public_kw_source()
  skip_if(is.null(path), "package sources unavailable")
  lines <- readLines(path, warn = FALSE)
  owner <- mpi_public_kw_source_region(
    lines,
    "^static int NP_OUTER_PACK_ADJACENT_HOT_ALIGN kernel_weighted_sum_np_ctx_ex\\(",
    "^int kernel_weighted_sum_np_ctx\\("
  )
  layout <- mpi_public_kw_source_region(
    lines,
    "^static int np_jksum_mpi_contiguous_row_layout\\(",
    "^#endif$"
  )
  compact <- gsub("[[:space:]]+", " ", owner)
  compact_layout <- gsub("[[:space:]]+", " ", layout)

  expect_match(
    compact,
    "kw_work = kw;",
    fixed = TRUE
  )
  expect_match(
    compact,
    paste0(
      "if(kw_equal_counts){ np_mpi_allgather_in_place_double( ",
      "igatherv[my_rank], kw_work, igatherv[0], ",
      "\"kernel weights MPI_Allgather\");"
    ),
    fixed = TRUE
  )
  expect_match(
    compact,
    paste0(
      "} else { np_mpi_allgatherv_in_place_double( ",
      "igatherv[my_rank], kw_work, igatherv, idisplsv, ",
      "\"kernel weights MPI_Allgatherv\");"
    ),
    fixed = TRUE
  )
  expect_match(
    compact_layout,
    "plane_elements > (size_t)INT_MAX",
    fixed = TRUE
  )
  expect_match(
    compact_layout,
    "rank_rows = MIN((size_t)row_stride, (size_t)total_rows - rank_start);",
    fixed = TRUE
  )

  expect_false(grepl("kw_work = (double *)calloc", owner, fixed = TRUE))
  expect_false(grepl("kw_work = (double *)malloc", owner, fixed = TRUE))
  expect_false(grepl("kw_recvcounts", owner, fixed = TRUE))
  expect_false(grepl("kw_displs", owner, fixed = TRUE))
  expect_false(grepl("memcpy(kw, kw_work", owner, fixed = TRUE))
  expect_false(grepl("free(kw_work)", owner, fixed = TRUE))

  expect_match(
    compact,
    paste0(
      "np_mpi_allreduce_in_place_double( ",
      "pkw_output->data, pkw_count, MPI_SUM, ",
      "\"derivative kernel weights MPI_Allreduce\");"
    ),
    fixed = TRUE
  )

  ## This owner must not bypass the progress-aware typed wrappers. Raw
  ## collectives remain valid only inside the wrapper fallback definitions.
  expect_false(grepl("MPI_Allgather[[:space:]]*\\(", owner, perl = TRUE))
  expect_false(grepl("MPI_Allgatherv[[:space:]]*\\(", owner, perl = TRUE))
  expect_false(grepl("MPI_Allreduce[[:space:]]*\\(", owner, perl = TRUE))
})
