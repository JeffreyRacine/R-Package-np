test_that("selector call equality ignores only the top-level source annotation", {
  same <- getFromNamespace(".np_bw_selector_call_identical", "npRmpi")
  plain <- quote(eval(call_obj, envir = value.env))
  annotated <- plain
  attr(annotated, "srcref") <- structure(1:8, class = "srcref")
  before <- serialize(annotated, NULL)
  expect_true(same(plain, plain))
  expect_true(same(plain, annotated))
  expect_true(same(annotated, plain))
  expect_identical(serialize(annotated, NULL), before)
  other <- annotated
  attr(other, "other") <- TRUE
  expect_false(same(plain, other))
  other <- annotated
  other[[2L]] <- quote(other_call)
  expect_false(same(plain, other))
  other <- annotated
  other[["extra"]] <- TRUE
  expect_false(same(plain, other))
  other <- annotated
  inner <- quote(f(call_obj))
  other[[2L]] <- inner
  a <- other
  attr(inner, "srcref") <- structure(1:8, class = "srcref")
  other[[2L]] <- inner
  expect_false(same(a, other))
  expect_false(same(plain, as.list(plain)))
})
