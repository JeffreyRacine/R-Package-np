test_that(".np_eval_call_arg prefers eval env and falls back to caller symbols", {
  call_obj <- as.call(list(as.name("f"), x = as.name("x")))
  eval_env <- new.env(parent = emptyenv())
  eval_env$x <- 11L
  environment(call_obj) <- eval_env

  caller_env <- new.env(parent = emptyenv())
  caller_env$x <- 22L
  expect_identical(.np_eval_call_arg(call_obj, "x", caller_env = caller_env), 11L)

  rm("x", envir = eval_env)
  expect_identical(.np_eval_call_arg(call_obj, "x", caller_env = caller_env), 22L)
})

test_that(".np_eval_bws_call_arg preserves bws fallback for missing symbols", {
  call_obj <- as.call(list(as.name("f"), txdat = as.name("txdat")))
  environment(call_obj) <- new.env(parent = emptyenv())

  bws <- list(call = call_obj, txdat = 1:4)
  expect_identical(.np_eval_bws_call_arg(bws, "txdat"), 1:4)
})
