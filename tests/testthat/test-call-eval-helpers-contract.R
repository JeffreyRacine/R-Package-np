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

test_that("call arg helpers preserve explicit NULL symbol values", {
  call_obj <- as.call(list(as.name("f"), x = as.name("x")))
  eval_env <- new.env(parent = emptyenv())
  eval_env$x <- NULL
  environment(call_obj) <- eval_env

  caller_env <- new.env(parent = emptyenv())
  caller_env$x <- 22L
  expect_null(.np_eval_call_arg(call_obj, "x", caller_env = caller_env))

  eval_env$txdat <- NULL
  call_obj2 <- as.call(list(as.name("f"), txdat = as.name("txdat")))
  environment(call_obj2) <- eval_env
  bws <- list(call = call_obj2, txdat = 1:4)
  expect_null(.np_eval_bws_call_arg(bws, "txdat"))
})

test_that(".np_eval_call_arg prefers explicit caller_env over ambient stack frames", {
  call_obj <- as.call(list(as.name("f"), x = as.name("x")))
  eval_env <- new.env(parent = emptyenv())
  environment(call_obj) <- eval_env

  caller_env <- new.env(parent = emptyenv())
  caller_env$x <- 22L

  out <- local({
    x <- 99L
    .np_eval_call_arg(call_obj, "x", caller_env = caller_env)
  })
  expect_identical(out, 22L)
})

test_that("call arg helpers route evaluation through shared frame helper", {
  call.body <- paste(deparse(body(.np_eval_call_arg), width.cutoff = 500L), collapse = " ")
  bws.body <- paste(deparse(body(.np_eval_bws_call_arg), width.cutoff = 500L), collapse = " ")

  expect_match(call.body, "\\.np_try_eval_in_frames\\(expr, eval_env = eval\\.env\\)")
  expect_match(bws.body, "\\.np_try_eval_in_frames\\(expr, eval_env = eval\\.env, search_frames = FALSE\\)")
})

test_that("explodePipe resolves formula symbols from provided env", {
  env <- new.env(parent = baseenv())
  env$fml <- y ~ x | z

  expect_identical(explodePipe(quote(fml), env = env),
                   explodePipe(env$fml, env = env))
})

test_that("explodePipe evaluates non-symbol formulas through strict frame helper", {
  fn.body <- paste(deparse(body(explodePipe), width.cutoff = 500L), collapse = " ")
  expect_match(fn.body, "\\.np_try_eval_in_frames\\(formula, eval_env = env, search_frames = FALSE\\)")
})
