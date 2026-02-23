test_that(".np_try_eval_in_frames resolves symbols from eval env first", {
  env <- new.env(parent = emptyenv())
  env$x <- 11L
  out <- .np_try_eval_in_frames(quote(x), eval_env = env)
  expect_true(out$ok)
  expect_identical(out$value, 11L)
})

test_that(".np_try_eval_in_frames resolves symbols from caller frames", {
  val <- local({
    x <- 17L
    .np_try_eval_in_frames(quote(x), eval_env = new.env(parent = emptyenv()))
  })
  expect_true(val$ok)
  expect_identical(val$value, 17L)
})
