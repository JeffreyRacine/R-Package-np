test_that("npscoefbw NOMAD service classifies malformed worker tasks", {
  worker <- getFromNamespace(".npscoefbw_nomad_slave_loop", "npRmpi")
  invalid_result <- getFromNamespace(".npscoefbw_nomad_invalid_worker_result", "npRmpi")
  task_error <- getFromNamespace(".npscoefbw_nomad_service_task_error", "npRmpi")
  worker_text <- paste(deparse(body(worker)), collapse = "\n")
  invalid_result_text <- paste(deparse(body(invalid_result)), collapse = "\n")
  task_error_text <- paste(deparse(body(task_error)), collapse = "\n")

  expect_true(grepl(".npscoefbw_nomad_invalid_worker_result", worker_text, fixed = TRUE))
  expect_true(grepl("malformed npscoefbw NOMAD service task", worker_text, fixed = TRUE))
  expect_true(grepl("unknown npscoefbw NOMAD service task", worker_text, fixed = TRUE))
  expect_true(grepl("malformed npscoefbw NOMAD eval task", worker_text, fixed = TRUE))
  expect_true(grepl("conditionMessage(e)", worker_text, fixed = TRUE))

  expect_true(grepl(".npscoefbw_nomad_service_task_error", invalid_result_text, fixed = TRUE))
  expect_true(grepl("invalid = 1L", invalid_result_text, fixed = TRUE))

  expect_true(grepl("npscoefbw.nomad.lp.service", task_error_text, fixed = TRUE))
  expect_true(grepl("task.error", task_error_text, fixed = TRUE))
  expect_true(grepl("task_type", task_error_text, fixed = TRUE))
})
