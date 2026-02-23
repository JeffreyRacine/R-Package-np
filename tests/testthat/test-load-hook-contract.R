test_that(".onLoad uses the shared dynamic-library helper", {
  load_body <- paste(deparse(body(.onLoad), width.cutoff = 500L), collapse = " ")
  expect_match(load_body, "\\.npRmpi_try_dynload\\(lib = lib, pkg = pkg\\)")
})

test_that(".npRmpi_try_dynload fails cleanly for unknown packages", {
  expect_false(.npRmpi_try_dynload(tempdir(), "npRmpi_definitely_missing_pkg"))
})
