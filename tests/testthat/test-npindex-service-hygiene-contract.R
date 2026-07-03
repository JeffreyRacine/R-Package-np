test_that("npindex service parameter predicate validates shape and finiteness", {
  param_ok <- getFromNamespace(".npindexbw_service_param_ok", "npRmpi")

  expect_true(param_ok(c(0.1, 0.2), 2L))
  expect_true(param_ok(c(0.1, 0.2, 0.3), 3L))
  expect_false(param_ok(c(0.1), 2L))
  expect_false(param_ok(c(0.1, Inf), 2L))
  expect_false(param_ok(c(0.1, NA_real_), 2L))
  expect_false(param_ok(c(0.1, 0.2), NA_integer_))
})

test_that("npindex service trace fields describe full-objective evaluation", {
  traced <- getFromNamespace(".npindexbw_eval_objective_service_traced", "npRmpi")
  body_text <- paste(deparse(body(traced)), collapse = "\n")

  expect_false(grepl("owned_rows", body_text, fixed = TRUE))
  expect_true(grepl("nominal_partition_rows", body_text, fixed = TRUE))
  expect_true(grepl("objective_rows", body_text, fixed = TRUE))
  expect_true(grepl("localize = FALSE", body_text, fixed = TRUE))
})

test_that("npindex Ichimura service malformed tasks resync via task errors", {
  worker <- getFromNamespace(".npindexbw_ichimura_lp_service_worker_loop", "npRmpi")
  body_text <- paste(deparse(body(worker)), collapse = "\n")

  expect_false(grepl("stop(\"malformed npindex Ichimura LP service task", body_text, fixed = TRUE))
  expect_false(grepl("stop(\"unknown npindex Ichimura LP service task", body_text, fixed = TRUE))
  expect_true(grepl(".npindexbw_ichimura_lp_service_task_error", body_text, fixed = TRUE))
  expect_true(grepl("unexpected npindex Ichimura LP service identifier", body_text, fixed = TRUE))
  expect_true(grepl("malformed npindex Ichimura LP eval task", body_text, fixed = TRUE))
  expect_true(grepl("unknown npindex Ichimura LP service task", body_text, fixed = TRUE))
})

test_that("npindex Klein-Spady service master and worker share eval predicate", {
  worker <- getFromNamespace(".npindexbw_kleinspady_lp_service_worker_loop", "npRmpi")
  master <- getFromNamespace(".npindexbw_kleinspady_lp_service_eval", "npRmpi")
  worker_text <- paste(deparse(body(worker)), collapse = "\n")
  master_text <- paste(deparse(body(master)), collapse = "\n")

  expect_true(grepl(".npindexbw_service_param_ok", worker_text, fixed = TRUE))
  expect_true(grepl("ncol(xmat)", worker_text, fixed = TRUE))
  expect_true(grepl(".npindexbw_service_param_ok", master_text, fixed = TRUE))
  expect_true(grepl("ncol(xmat)", master_text, fixed = TRUE))
})
