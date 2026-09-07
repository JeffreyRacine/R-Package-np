test_that("restoration notices use the existing owner and scheduled render", {
  # A private lexical copy keeps the namespace/global progress owner untouched.
  scope <- new.env(parent = environment(.np_progress_begin))
  for (name in ls(environment(.np_progress_begin), all.names = TRUE)) {
    if (!startsWith(name, ".np_progress_")) next
    value <- get(name, environment(.np_progress_begin))
    if (is.function(value)) environment(value) <- scope
    assign(name, value, scope)
  }
  scope$.np_progress_runtime <- new.env(parent = emptyenv())
  scope$.np_progress_registry <- scope$.np_progress_make_registry()
  scope$.np_progress_is_interactive <- function() TRUE
  scope$.np_progress_bandwidth_worker_silent <- function() FALSE
  scope$.np_progress_renderer_for_surface <- function(...) "single_line"
  scope$.np_progress_resolve_message_muffling <- identity
  clock <- 0; width <- 80L; probes <- 0L; events <- list()
  scope$.np_progress_now <- function() clock
  scope$.np_progress_output_width <- function() {probes <<- probes+1L; width}
  scope$.np_progress_render_single_line <- function(snapshot, event) {
    events[[length(events)+1L]] <<- list(snapshot=snapshot, event=event)
  }
  for (enhanced in c(FALSE, TRUE)) for (forwarded in c(FALSE, TRUE)) {
    scope$.np_progress_reset_registry(); events <- list(); clock <- 0
    state <- scope$.np_progress_begin("Bandwidth selection", domain="bandwidth")
    state$enabled <- state$visible <- TRUE
    state$start_note_pending <- FALSE
    state$last_emit <- 0; state$throttle_sec <- 1
    state$last_done <- 17L
    state$bandwidth_progress_common <- enhanced
    state$bandwidth_mode <- "iteration"
    slot <- if (forwarded) "bandwidth_forward_state" else "bandwidth_state"
    scope$.np_progress_runtime$bandwidth_forward_active <- forwarded
    scope$.np_progress_runtime[[slot]] <- state
    labels <- c("start 2 retry 0: valid start restored", "start 2 restored", "start restored")
    before <- probes
    scope$.np_progress_bandwidth_notice(labels)
    expect_identical(probes, before)
    expect_length(events, 0L)
    clock <- .5; scope$.np_progress_bandwidth_activity_step()
    expect_length(events, 0L)
    expect_identical(scope$.np_progress_runtime[[slot]]$bandwidth_notice, labels)
    clock <- 1.5; scope$.np_progress_bandwidth_activity_step()
    expect_length(events, 1L)
    expect_match(events[[1L]]$snapshot$render_line, "restored")
    expect_identical(events[[1L]]$snapshot$id, state$id)
    expect_identical(scope$.np_progress_runtime[[slot]]$last_done, 17L)
    expect_null(scope$.np_progress_runtime[[slot]]$bandwidth_notice)
    expect_identical(probes, before+1L)
    clock <- 3; scope$.np_progress_bandwidth_activity_step()
    expect_false(grepl("restored", events[[2L]]$snapshot$render_line))
    scope$.np_progress_end(scope$.np_progress_runtime[[slot]])
    expect_identical(tail(events,1L)[[1L]]$event, "finish")
    expect_null(scope$.np_progress_registry$active_id)
  }
  scope$.np_progress_runtime$bandwidth_forward_active <- FALSE
  for (kind in c("disabled", "invisible", "worker", "absent")) {
    state$enabled <- kind != "disabled"; state$visible <- kind != "invisible"
    scope$.np_progress_runtime$bandwidth_state <- if(kind=="absent") NULL else state
    scope$.np_progress_bandwidth_worker_silent <- function() kind=="worker"
    scope$.np_progress_bandwidth_notice(labels)
    expect_null(scope$.np_progress_runtime$bandwidth_state$bandwidth_notice)
  }
})

test_that("restoration wording fits the owner's existing width budget", {
  for (prefix in c("[np]", "[npRmpi]")) for (width in 20:120) {
    text <- .np_progress_bandwidth_notice_line(
      paste(prefix,"Bandwidth selection (iteration 100000, elapsed 9999.0s)"),
      prefix, c("start 123456 retry 100: valid start restored", "start 123456 restored", "start restored"),
      width)
    expect_true(nchar(text, type="width") <= width)
    expect_match(text, "restored")
  }
})
