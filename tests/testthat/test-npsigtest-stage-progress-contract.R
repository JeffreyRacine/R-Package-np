test_that("streamed native tiles forward activity and unwind progress safely", {
  
  ns <- asNamespace("np")
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(53)
  x <- data.frame(x = rnorm(40), f = factor(rep(letters[1:2], 20)))
  y <- x$x + rnorm(40)
  bw <- npregbw(xdat = x, ydat = y, bws = c(.6, .2),
                bandwidth.compute = FALSE)
  tile <- get(".np_npsig_streamed_iid_tile", ns)
  responses <- cbind(y, rev(y))
  run <- function(index = 1L, pivotal = TRUE) tile(bw, x, index,
    response.matrix = responses, null.mean = y, residual.pool = y,
    pivotal = pivotal)
  reference <- run()
  pulses <- 0L
  prior.forward <- get(".np_progress_runtime", ns)$fit_forward
  bindings <- list(
    .np_progress_is_interactive = function() TRUE,
    .np_fit_progress_step = function(...) {
      pulses <<- pulses + 1L
      invisible(NULL)
    }
  )
  with_np_progress_bindings(bindings, {
    options(np.messages = TRUE)
    expect_identical(run(), reference)
    expect_gt(pulses, 0L)
    pulses <- 0L
    invisible(run(1L, FALSE))
    expect_gt(pulses, 0L)
    pulses <- 0L
    invisible(run(2L, FALSE))
    expect_gt(pulses, 0L)
    expect_error(tile(bw, x, 1L, response.matrix = responses * Inf,
      null.mean = y, residual.pool = y, pivotal = TRUE), "failed")
    expect_identical(run(), reference)
    expect_null(get(".np_progress_runtime", ns)$fit_state)
    expect_identical(get(".np_progress_runtime", ns)$fit_forward, prior.forward)
    options(np.messages = FALSE)
    pulses <- 0L
    expect_identical(run(), reference)
    expect_identical(pulses, 0L)
  })
})

test_that("streamed individual preparation does not refit unrestricted residuals", {
  
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(72)
  x <- data.frame(x = rnorm(45), z = factor(rep(1:3, 15)))
  y <- x$x + rnorm(45)
  bw <- npregbw(xdat = x, ydat = y, bws = c(.7, .25),
                bandwidth.compute = FALSE)
  original <- getFromNamespace(".npreg_complete", "np")
  calls <- list()
  with_np_progress_bindings(list(.npreg_complete = function(...) {
    args <- list(...)
    calls[[length(calls) + 1L]] <<- args
    original(...)
  }), {
    result <- npsigtest(bw, xdat = x, ydat = y, B = 9, index = 1L)
    expect_s3_class(result, "sigtest")
    joint <- npsigtest(bw, xdat = x, ydat = y, B = 9, joint = TRUE)
    expect_s3_class(joint, "sigtest")
  })
  expect_false(any(vapply(calls, function(args) isTRUE(args$residuals), logical(1))))
  expect_true(any(vapply(calls, function(args) isTRUE(args$gradients), logical(1))))
})
