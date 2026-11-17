test_that("untagging untagged inputs preserves vector sharing through serialization", {
  skip_if_not(capabilities("profmem"))
  d <- data.frame(x = sin(seq_len(128L)), y = cos(seq_len(128L)))
  before <- vapply(d, tracemem, character(1L))
  on.exit(lapply(d, untracemem), add = TRUE)
  out <- .npRmpi_autodispatch_untag(list(training = d, call = quote(f(x))))
  invisible(serialize(out, NULL))
  after <- vapply(out$training, tracemem, character(1L))
  on.exit(lapply(out$training, untracemem), add = TRUE)
  expect_identical(after, before)
  expect_identical(out$training, d)
})

test_that("untagging removes real dispatch tags and preserves unrelated structure", {
  tags <- c("npRmpi.dispatch.mode", "npRmpi.autodispatch.remote",
            "npRmpi.autodispatch.fingerprint", "npRmpi.autodispatch.lease")
  d <- data.frame(x = 1:8, g = factor(rep(c("a", "b"), 4)))
  attr(d$x, "label") <- "keep"
  clean <- list(training = d, empty = NULL, call = quote(f(, x, NULL)), formula = y ~ x)
  tagged <- clean
  for (name in tags) {
    attr(tagged, name) <- "remove"
    attr(tagged$training$x, name) <- "remove"
  }
  expect_identical(.npRmpi_autodispatch_untag(tagged), clean)
  expect_identical(attr(tagged$training$x, "label", exact = TRUE), "keep")
  for (name in tags) expect_identical(attr(tagged, name, exact = TRUE), "remove")
})
