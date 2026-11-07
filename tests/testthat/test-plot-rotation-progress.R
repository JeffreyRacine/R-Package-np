test_that("rotation measures remaining frames independently of inference", {
  old <- options(np.messages=TRUE, np.plot.progress=TRUE)
  on.exit(options(old), add=TRUE)
  now <- 0
  z <- capture_progress_shadow_trace(np:::.np_plot_progress_run({
    np:::.np_plot_progress_plan(1L)
    target <- np:::.np_plot_progress_target_begin()
    now <- 26
    np:::.np_plot_progress_notify(done=9999L,total=9999L,force=TRUE)
    np:::.np_plot_progress_target_end(target)
    previous <- target$phase
    now <- 30
    rotation <- np:::.np_plot_rotation_progress_begin(64L)
    now <- 32
    rotation <- np:::.np_plot_rotation_progress_tick(rotation,16L)
    now <- 34
    rotation <- np:::.np_plot_rotation_progress_tick(rotation,32L)
    expect_identical(target$completed,1L)
    expect_identical(target$total,1L)
    np:::.np_plot_rotation_progress_end(rotation)
    expect_identical(target$phase,previous)
    expect_null(target$rep)
  },enabled=TRUE),force_renderer="single_line",now=function()now)
  lines <- vapply(z$trace,`[[`,"","line")
  expect_true(any(grepl("frames 16/64.*elapsed 32.0s, frame eta 6.0s.*Rotating plot",lines)))
  expect_true(any(grepl("frames 32/64.*elapsed 34.0s, frame eta 4.0s",lines)))
  expect_false(any(grepl("frames (16|32)/64.*eta finishing",lines)))
  expect_length(unique(vapply(z$trace,`[[`,"","id")),1L)
  expect_null(np:::.np_plot_progress_runtime$context)
})

test_that("rotation restores its parent phase on error without success", {
  old <- options(np.messages=TRUE,np.plot.progress=TRUE)
  on.exit(options(old),add=TRUE)
  z <- capture_progress_shadow_trace({
    expect_error(np:::.np_plot_progress_run({
      np:::.np_plot_progress_outer_run({
        context <- np:::.np_plot_progress_runtime$context
        previous <- context$phase
        rotate <- function() {
          state <- np:::.np_plot_rotation_progress_begin(64L)
          on.exit(np:::.np_plot_rotation_progress_end(state),add=TRUE)
          stop("renderer sentinel")
        }
        expect_error(rotate(),"renderer sentinel")
        expect_identical(context$phase,previous)
        stop("outer sentinel")
      })
    },enabled=TRUE),"outer sentinel")
  },force_renderer="single_line")
  expect_false(any(vapply(z$trace,function(x)x$event=="finish",TRUE)))
  expect_null(np:::.np_plot_progress_runtime$context)
  expect_null(np:::.np_progress_registry$active_id)
})
