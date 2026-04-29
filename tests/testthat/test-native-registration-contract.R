test_that("native registration disables dynamic lookup and registers MPI entry points", {
  src <- test_path("..", "..", "src", "np_init.c")
  skip_if_not(file.exists(src), "source C files unavailable in installed test context")
  lines <- readLines(src, warn = FALSE)
  text <- paste(lines, collapse = "\n")

  expect_match(text, "R_useDynamicSymbols\\(dll, FALSE\\);")
  expect_match(text, "\\{\"mpi_finalize\",\\s*\\(DL_FUNC\\) &mpi_finalize,\\s*0\\}")
  expect_match(text, "\\{\"mpi_gather\",\\s*\\(DL_FUNC\\) &mpi_gather,\\s*5\\}")
  expect_match(text, "\\{\"mpi_sendrecv\",\\s*\\(DL_FUNC\\) &mpi_sendrecv,\\s*10\\}")
  expect_match(text, "\\{\"mpi_wait\",\\s*\\(DL_FUNC\\) &mpi_wait,\\s*2\\}")
})

test_that("package metadata declares the MPI system dependency", {
  desc <- test_path("..", "..", "DESCRIPTION")
  skip_if_not(file.exists(desc), "DESCRIPTION unavailable in installed test context")
  lines <- readLines(desc, warn = FALSE)
  expect_true(any(grepl("^SystemRequirements:\\s*MPI\\s*$", lines)))
})
