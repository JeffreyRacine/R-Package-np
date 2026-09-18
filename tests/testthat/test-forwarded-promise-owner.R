test_that("only an exact activation's formals own matched forwarded arguments", {
  resolve <- getFromNamespace(".npRmpi_autodispatch_resolve_owned_arg", "npRmpi")
  dots.only <- function(...) {
    ydat <- rep(11, 3)
    resolve(quote(..1), "ydat", environment(), "")
  }
  formal <- function(ydat, ...) {
    # ..1 has been consumed by S3 matching; the original formal now owns it.
    resolve(quote(..1), "ydat", environment(), "")
  }
  named <- function(...) {
    ydat <- 11
    resolve(quote(..1), "ydat", environment(), "ydat")
  }
  expect_identical(dots.only(1:3), 1:3)
  expect_identical(formal(1:3, 88), 1:3)
  expect_identical(named(ydat = 1:3), 1:3)
  n <- 0L
  expect_identical(dots.only({ n <- n + 1L; 1:3 }), 1:3)
  expect_identical(n, 1L)
  expect_error(dots.only({ n <- n + 1L; stop("original promise") }), "original promise")
  expect_identical(n, 2L)
  targets <- getFromNamespace(".npRmpi_autodispatch_target_args", "npRmpi")()
  sweep <- function(name, ...) {
    assign(name, "unrelated local", envir = environment())
    resolve(quote(..1), name, environment(), "")
  }
  # name is a control formal in this test, not a materialized package argument.
  for (target in setdiff(targets, "name"))
    expect_identical(sweep(target, 1:3), 1:3, info = target)
})
