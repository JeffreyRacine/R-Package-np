test_that("plot heartbeat retains completed work between render checkpoints", {
  ns <- asNamespace("npRmpi")
  runtime <- get(".np_plot_progress_runtime", ns)
  previous <- runtime$context
  on.exit(runtime$context <- previous)
  context <- new.env(parent = emptyenv())
  context$rep <- NULL
  context$stage <- "preparing"
  context$fraction <- 0
  runtime$context <- context
  record <- get(".np_plot_progress_record", ns)
  record(done=8L,total=399L)
  expect_identical(context$rep,"rep 8/399")
  expect_equal(context$fraction,8/399)
  record(stage="preparing operator",done=4L,total=40L,unit="row")
  expect_identical(context$rep,"row 4/40")
  expect_equal(context$fraction,8/399)
  record(done=2L,total=4L,unit="block")
  expect_identical(context$rep,"block 2/4")
  expect_equal(context$fraction,.5)
  record()
  expect_identical(context$rep,"block 2/4")
  # Exercise the early-return path: no visible progress or clock decision
  # should be needed to retain newly completed work for the outer heartbeat.
  tick <- get(".np_plot_progress_tick", ns)
  env <- new.env(parent=ns)
  env$.np_progress_now <- function() 1
  env$.np_progress_maybe_emit_start_note <- function(state,now) state
  environment(tick) <- env
  state <- list(total=399L,plot_context=context,plot_unit="rep",
                checkpoints=399L,next_checkpoint_idx=1L,
                last_emitted_done=0L,start_note_pending=FALSE,
                last_emit=1,throttle_sec=100)
  result <- tick(state,16L)
  expect_identical(result$last_done,16L)
  expect_identical(context$rep,"rep 16/399")
  expect_identical(result$last_emit,1)
})
