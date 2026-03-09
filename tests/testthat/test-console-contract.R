new_line_console <- getFromNamespace("newLineConsole", "np")
to_msg <- getFromNamespace("toMsg", "np")

test_that("toMsg tab expansion remains deterministic and side-effect free", {
  console <- new_line_console()
  out <- to_msg("a\tb", console = console)

  expect_identical(out$msg, paste0("a", strrep(" ", 6L), "b", strrep(" ", 7L)))
  expect_identical(out$len, nchar(out$msg))
  expect_false(grepl("<<-",
                     paste(deparse(body(to_msg), width.cutoff = 500L), collapse = " "),
                     fixed = TRUE))
})
