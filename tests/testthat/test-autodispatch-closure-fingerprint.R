# Focused, non-MPI contract: ordinary R transport changes the environment field
# of evaluated promises captured by closures. Keep this separate from lifecycle
# matrices; the original public saved-formula test protects integration.
test_that("local closure fingerprints survive promise transport without losing captures", {
  fingerprint <- getFromNamespace(".npRmpi_autodispatch_fingerprint", "npRmpi")
  make <- eval(quote(function(value) {
    force(value)
    function(x) x + value
  }), baseenv())
  f <- make(2)
  x <- list(callback = f)
  owner <- environment(f)
  before <- serialize(x, NULL)
  expect_identical(fingerprint(x), fingerprint(unserialize(before)))
  expect_identical(serialize(x, NULL), before)
  expect_identical(environment(f), owner)
  expect_false(identical(fingerprint(x), fingerprint(list(callback = make(3)))))
  old <- fingerprint(x)
  assign("value", 3, envir = owner)
  expect_false(identical(old, fingerprint(x)))
  expect_identical(f(1), 4)
  attr(x, "callback") <- make(4)
  expect_identical(fingerprint(x), fingerprint(unserialize(serialize(x, NULL))))
})

test_that("hash snapshots do not evaluate lazy or active captured bindings", {
  fingerprint <- getFromNamespace(".npRmpi_autodispatch_fingerprint", "npRmpi")
  make <- eval(quote(function(value) function(x) x + value), baseenv())
  lazy <- make(stop("hash forced an argument"))
  expect_identical(fingerprint(list(callback = lazy)),
                   fingerprint(unserialize(serialize(list(callback = lazy), NULL))))
  expect_error(lazy(1), "hash forced an argument")
  owner <- new.env(parent = baseenv())
  makeActiveBinding("value", function(v) stop("hash read an active binding"), owner)
  f <- eval(quote(function(x) x + value), owner)
  expect_identical(fingerprint(list(callback = f)),
                   fingerprint(unserialize(serialize(list(callback = f), NULL))))
  expect_true(bindingIsActive("value", owner))
  expect_identical(environment(f), owner)
  expect_error(f(1), "hash read an active binding")
})

test_that("ordinary and context fingerprint representations stay unchanged", {
  fingerprint <- getFromNamespace(".npRmpi_autodispatch_fingerprint", "npRmpi")
  digest <- getFromNamespace(".npRmpi_autodispatch_raw_md5", "npRmpi")
  expected <- function(x) {
    bytes <- serialize(x, NULL, xdr = FALSE)
    sprintf("md5:%d:%s", length(bytes), digest(bytes))
  }
  for (x in list(list(a = 1:10), list(callback = stats::na.exclude),
                 list(callback = "na.exclude"), list(callback = NULL)))
    expect_identical(fingerprint(x), expected(x))
  make <- eval(quote(function(value) { force(value); function(x) x + value }), baseenv())
  f <- make(2)
  # Source references have their own incumbent attribute normalization.
  attributes(f) <- NULL
  x <- list(callback = f)
  expect_identical(fingerprint(x, "context.v1"), expected(x))
})
