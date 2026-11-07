test_that("manual RNG refusal handling is confined to receiving commands", {
  context <- getFromNamespace(".npRmpi_with_manual_bcast_context", "npRmpi")
  refusal <- structure(list(message = "agreed test refusal", call = NULL),
                       class = c("np_rng_scope_refusal", "error", "condition"))
  before <- options()[c("npRmpi.manual.bcast.context", "npRmpi.autodispatch.disable")]
  on.exit(options(before), add = TRUE)
  for (receive in c(FALSE, TRUE)) {
    for (value in list(NULL, 1L, simpleError("returned"), structure("x", class="try-error"))) {
      expect_identical(withVisible(context(value, receive)), withVisible(value))
      expect_identical(withVisible(context(invisible(value), receive)), withVisible(invisible(value)))
    }
    expect_error(context(stop("ordinary"), receive), "ordinary")
    expect_warning(context(warning("ordinary warning"), receive), "ordinary warning")
    interrupt <- structure(list(message="test", call=NULL), class=c("interrupt", "condition"))
    expect_identical(tryCatch(context(stop(interrupt), receive), interrupt=identity), interrupt)
  }
  cleanup <- FALSE; after <- FALSE
  body <- function() { on.exit(cleanup <<- TRUE); stop(refusal) }
  expect_identical(context({context(body()); after <- TRUE}, TRUE), refusal)
  expect_true(cleanup); expect_false(after)
  expect_identical(tryCatch(context(stop(refusal)), error=identity), refusal)
  expect_identical(options()[names(before)], before)
})
