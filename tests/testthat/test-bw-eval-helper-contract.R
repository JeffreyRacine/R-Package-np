np_try_eval_in_frames <- getFromNamespace(".np_try_eval_in_frames", "np")
np_eval_bw_call <- getFromNamespace(".np_eval_bw_call", "np")

test_that(".np_try_eval_in_frames resolves symbols from eval env first", {
  env <- new.env(parent = emptyenv())
  env$x <- 11L
  out <- np_try_eval_in_frames(quote(x), eval_env = env)
  expect_true(out$ok)
  expect_identical(out$value, 11L)
})

test_that(".np_try_eval_in_frames preserves NULL symbol values", {
  env <- new.env(parent = emptyenv())
  env$x <- NULL
  out <- np_try_eval_in_frames(quote(x), eval_env = env)
  expect_true(out$ok)
  expect_null(out$value)
})

test_that(".np_try_eval_in_frames resolves symbols from caller frames", {
  val <- local({
    x <- 17L
    np_try_eval_in_frames(quote(x), eval_env = new.env(parent = emptyenv()))
  })
  expect_true(val$ok)
  expect_identical(val$value, 17L)
})

test_that(".np_try_eval_in_frames can disable caller frame fallback", {
  out <- local({
    x <- 19L
    np_try_eval_in_frames(quote(x), eval_env = new.env(parent = emptyenv()), search_frames = FALSE)
  })
  expect_false(out$ok)
  expect_true(inherits(out$error, "error"))
})

test_that(".np_try_eval_in_frames returns an error object when resolution fails", {
  out <- np_try_eval_in_frames(quote(np_missing_symbol_contract), eval_env = new.env(parent = emptyenv()))
  expect_false(out$ok)
  expect_true(inherits(out$error, "error"))
  expect_match(conditionMessage(out$error), "np_missing_symbol_contract")
})

test_that(".np_try_eval_in_frames evaluates non-symbol expressions in caller frames", {
  out <- local({
    x <- 4L
    np_try_eval_in_frames(quote(x + 1L), eval_env = new.env(parent = emptyenv()))
  })
  expect_true(out$ok)
  expect_identical(out$value, 5L)
})

test_that(".np_try_eval_in_frames can disable caller fallback for non-symbol expressions", {
  out <- local({
    x <- 4L
    np_try_eval_in_frames(
      quote(x + 1L),
      eval_env = new.env(parent = emptyenv()),
      search_frames = FALSE
    )
  })
  expect_false(out$ok)
  expect_true(inherits(out$error, "error"))
})

test_that(".np_try_eval_in_frames honors enclos for non-symbol expressions", {
  eval_env <- list()
  enclos <- new.env(parent = baseenv())
  enclos$z <- 8L
  out <- np_try_eval_in_frames(
    quote(z + 2L),
    eval_env = eval_env,
    enclos = enclos,
    search_frames = FALSE
  )
  expect_true(out$ok)
  expect_identical(out$value, 10L)
})

test_that(".np_eval_bw_call reports underlying evaluation errors", {
  call_obj <- as.call(list(as.name("identity"), as.name("np_missing_bw_arg_contract")))
  environment(call_obj) <- new.env(parent = baseenv())

  expect_error(
    np_eval_bw_call(call_obj, caller_env = new.env(parent = baseenv())),
    "np_missing_bw_arg_contract"
  )
})
