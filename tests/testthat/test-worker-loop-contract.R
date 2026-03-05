test_that("attach worker loop delegates to shared worker loop helper", {
  fn <- getFromNamespace(".npRmpi_session_attach_worker_loop", "npRmpi")
  fn.body <- paste(deparse(body(fn), width.cutoff = 500L), collapse = " ")
  expect_match(fn.body, "\\.npRmpi_worker_loop\\(")
})

test_that("spawn/profile worker scripts use shared worker loop helper", {
  spawn.path <- system.file("slavedaemon.R", package = "npRmpi")
  expect_true(nzchar(spawn.path) && file.exists(spawn.path))
  spawn.lines <- readLines(spawn.path, warn = FALSE)
  expect_true(any(grepl("\\.npRmpi_worker_loop\\(", spawn.lines)))

  profile.path <- system.file("Rprofile", package = "npRmpi")
  expect_true(nzchar(profile.path) && file.exists(profile.path))
  profile.lines <- readLines(profile.path, warn = FALSE)
  expect_true(any(grepl("\\.npRmpi_worker_loop\\(", profile.lines)))
})

test_that("worker message normalizer canonicalizes do.call payloads", {
  norm <- getFromNamespace(".npRmpi_worker_normalize_message", "npRmpi")
  env <- new.env(parent = emptyenv())
  msg <- as.call(list(as.name("do.call"), as.name("mean"), quote(list(1, 2, 3)), env))
  out <- norm(msg)
  expect_true(is.call(out))
  expect_identical(as.character(out[[1L]]), "do.call")
  expect_identical(out[[4L]], FALSE)
  expect_true(is.environment(out[[5L]]))
})

test_that("worker break detector handles raw and wrapped shutdown sentinels", {
  is_break <- getFromNamespace(".npRmpi_worker_is_break_message", "npRmpi")

  expect_true(is_break("kaerb"))
  wrapped <- as.call(list(as.name(".npRmpi_with_manual_bcast_context"), "kaerb"))
  expect_true(is_break(wrapped))
  wrapped_fun <- as.call(list(function(x) x, "kaerb"))
  expect_true(is_break(wrapped_fun))

  expect_false(is_break("not_break"))
  expect_false(is_break(quote(.npRmpi_with_manual_bcast_context(1 + 1))))
})
