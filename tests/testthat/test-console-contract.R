test_that("toMsg tab expansion remains deterministic and side-effect free", {
  console <- newLineConsole()
  out <- toMsg("a\tb", console = console)

  expect_identical(out$msg, paste0("a", strrep(" ", 6L), "b", strrep(" ", 7L)))
  expect_identical(out$len, nchar(out$msg))
  expect_false(grepl("<<-",
                     paste(deparse(body(toMsg), width.cutoff = 500L), collapse = " "),
                     fixed = TRUE))
})
