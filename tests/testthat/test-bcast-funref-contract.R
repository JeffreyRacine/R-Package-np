test_that(".npRmpi_bcast_cmd_funref resolves standard command heads", {
  expect_identical(.npRmpi_bcast_cmd_funref(quote(assign)), "assign")
  expect_identical(.npRmpi_bcast_cmd_funref("assign"), "assign")
  expect_identical(.npRmpi_bcast_cmd_funref(base::assign), base::assign)
})

test_that(".npRmpi_bcast_cmd_funref resolves namespace-qualified heads", {
  ref <- .npRmpi_bcast_cmd_funref(quote(base::assign))
  expect_true(is.function(ref))
  expect_identical(ref, base::assign)

  ref_get <- .npRmpi_bcast_cmd_funref(quote(base::get))
  expect_true(is.function(ref_get))
  expect_identical(ref_get, base::get)

  ref_internal <- .npRmpi_bcast_cmd_funref(quote(base:::print.default))
  expect_true(is.function(ref_internal))
  expect_identical(ref_internal, get("print.default", envir = asNamespace("base"), mode = "function", inherits = FALSE))

  ref_call_head <- .npRmpi_bcast_cmd_funref(quote((base::assign)("tmp", 1L)))
  expect_true(is.function(ref_call_head))
  expect_identical(ref_call_head, base::assign)
})

test_that(".npRmpi_bcast_cmd_funref no longer evals full command expressions", {
  fn.body <- paste(deparse(body(.npRmpi_bcast_cmd_funref), width.cutoff = 500L), collapse = " ")
  expect_false(grepl("eval\\(scmd", fn.body))
  expect_false(grepl("eval\\(hd", fn.body))
  expect_match(fn.body, "\\.npRmpi_eval_scmd\\(hd, envir = eval_env\\)")
})
