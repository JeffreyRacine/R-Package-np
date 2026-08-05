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

test_that("MPI public kernel weights collect into one native plane", {
  path <- locate_mpi_public_kw_source()
  skip_if(is.null(path), "package sources unavailable")
  source <- paste(readLines(path, warn = FALSE), collapse = "\n")
  compact <- gsub("[[:space:]]+", " ", source)

  expect_match(
    compact,
    "kw_work = kw;",
    fixed = TRUE
  )
  expect_match(
    compact,
    paste0(
      "if(kw_equal_counts){ MPI_Allgather(MPI_IN_PLACE, ",
      "igatherv[my_rank], MPI_DOUBLE, kw_work, igatherv[0], ",
      "MPI_DOUBLE, comm[1]);"
    ),
    fixed = TRUE
  )
  expect_match(
    compact,
    paste0(
      "} else { MPI_Allgatherv(MPI_IN_PLACE, ",
      "igatherv[my_rank], MPI_DOUBLE, kw_work, igatherv, idisplsv, ",
      "MPI_DOUBLE, comm[1]);"
    ),
    fixed = TRUE
  )
  expect_match(
    compact,
    "plane_elements > (size_t)INT_MAX",
    fixed = TRUE
  )
  expect_match(
    compact,
    "rank_rows = MIN((size_t)row_stride, (size_t)total_rows - rank_start);",
    fixed = TRUE
  )

  expect_false(grepl("kw_work = (double *)calloc", source, fixed = TRUE))
  expect_false(grepl("kw_work = (double *)malloc", source, fixed = TRUE))
  expect_false(grepl("kw_recvcounts", source, fixed = TRUE))
  expect_false(grepl("kw_displs", source, fixed = TRUE))
  expect_false(grepl("memcpy(kw, kw_work", source, fixed = TRUE))
  expect_false(grepl("free(kw_work)", source, fixed = TRUE))

  expect_match(
    compact,
    paste0(
      "MPI_Allreduce(MPI_IN_PLACE, pkw_output->data, pkw_count, ",
      "MPI_DOUBLE, MPI_SUM, comm[1]);"
    ),
    fixed = TRUE
  )
})
