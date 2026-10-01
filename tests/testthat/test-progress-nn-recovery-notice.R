nn_notice_scope <- function() {
  parent <- environment(.np_progress_begin)
  scope <- new.env(parent = parent)
  for (name in ls(parent, all.names = TRUE)) {
    if (!startsWith(name, ".np_progress_") && !startsWith(name, ".np_nn_")) next
    value <- get(name, parent)
    if (is.function(value)) environment(value) <- scope
    assign(name, value, scope)
  }
  scope$.np_progress_runtime <- new.env(parent = emptyenv())
  scope$.np_progress_registry <- scope$.np_progress_make_registry()
  scope$.np_progress_is_interactive <- function() TRUE
  scope$.np_progress_bandwidth_worker_silent <- function() FALSE
  scope$.np_progress_renderer_for_surface <- function(...) "single_line"
  scope$.np_progress_resolve_message_muffling <- identity
  scope$.np_progress_output_width <- function() 180L
  scope
}

test_that("NN recovery is an immediate one-time notice without search state changes", {
  s <- nn_notice_scope()
  events <- list(); clock <- 0
  s$.np_progress_now <- function() clock
  s$.np_progress_render_single_line <- function(snapshot, event)
    events[[length(events)+1L]] <<- list(snapshot=snapshot, event=event)
  for (forwarded in c(FALSE, TRUE)) for (total in c(1L, 10L)) {
    s$.np_progress_reset_registry(); events <- list(); clock <- 0
    state <- s$.np_progress_bandwidth_initialize_state(
      s$.np_progress_begin("Bandwidth selection"))
    state$enabled <- state$visible <- TRUE
    state$bandwidth_nmulti_total <- total
    state$bandwidth_multistart_current <- total
    state$bandwidth_multistart_completed <- total
    state$last_done <- 64L
    slot <- if (forwarded) "bandwidth_forward_state" else "bandwidth_state"
    s$.np_progress_runtime$bandwidth_forward_active <- forwarded
    s$.np_progress_runtime[[slot]] <- state
    set.seed(211); rng <- .Random.seed
    expect_true(s$.np_progress_signal_from_c("bandwidth_nn_recovery", "bandwidth"))
    expect_identical(.Random.seed, rng)
    expect_length(events, 1L)
    expect_match(events[[1L]]$snapshot$render_line,
                 "NN bandwidth infeasible; attempting recovery", fixed=TRUE)
    expect_false(grepl("100[.]0%|eta 0[.]0s|restored", events[[1L]]$snapshot$render_line))
    after <- s$.np_progress_runtime[[slot]]
    fields <- c("id","started","last_done","bandwidth_nmulti_total",
                "bandwidth_multistart_current","bandwidth_multistart_completed",
                "progress_provider","unknown_total_fields")
    expect_identical(after[fields], state[fields])
    expect_null(after$bandwidth_notice)
    clock <- 3
    s$.np_progress_bandwidth_activity_step(force=TRUE)
    expect_false(grepl("infeasible|recovery|recovering", tail(events,1L)[[1L]]$snapshot$render_line))
    s$.np_progress_end(s$.np_progress_runtime[[slot]])
    expect_identical(tail(events,1L)[[1L]]$event, "finish")
    expect_null(s$.np_progress_registry$active_id)
  }
})

test_that("NN recovery notice is silent without a visible nonworker owner", {
  s <- nn_notice_scope(); events <- list()
  s$.np_progress_render_single_line <- function(snapshot,event)
    events[[length(events)+1L]] <<- snapshot
  for (kind in c("disabled","invisible","worker","absent","not-owner")) {
    s$.np_progress_reset_registry()
    state <- s$.np_progress_begin("Bandwidth selection")
    state$enabled <- kind != "disabled"; state$visible <- kind != "invisible"
    s$.np_progress_runtime$bandwidth_state <- if(kind=="absent") NULL else state
    s$.np_progress_bandwidth_worker_silent <- function() kind=="worker"
    if(kind=="not-owner") s$.np_progress_registry$active_id <- "another-owner"
    before <- s$.np_progress_runtime$bandwidth_state
    s$.np_progress_bandwidth_nn_recovery()
    expect_identical(s$.np_progress_runtime$bandwidth_state, before)
  }
  expect_length(events, 0L)
})

test_that("nested recovery retains the canonical owner and its refinement provider", {
  s <- nn_notice_scope(); events <- list()
  s$.np_progress_render_single_line <- function(snapshot,event)
    events[[length(events)+1L]] <<- list(snapshot=snapshot,event=event)
  state <- s$.np_progress_begin("Refining bandwidth", domain="bandwidth")
  state$progress_provider <- "powell_refine"
  state$unknown_total_fields <- function(...) "refinement fields"
  state$bandwidth_notice <- c("ordinary restoration", "restoration")
  s$.np_progress_runtime$bandwidth_state <- state
  next.id <- s$.np_progress_registry$next_id
  s$.np_progress_with_nested_bandwidth_heartbeat({
    s$.np_progress_with_nested_bandwidth_heartbeat({
      s$.np_progress_bandwidth_nn_recovery()
    })
  })
  after <- s$.np_progress_runtime$bandwidth_state
  expect_identical(s$.np_progress_registry$next_id,next.id)
  expect_identical(s$.np_progress_registry$active_id,state$id)
  expect_identical(after[c("id","label","progress_provider","unknown_total_fields","started",
                          "bandwidth_notice")],
                   state[c("id","label","progress_provider","unknown_total_fields","started",
                           "bandwidth_notice")])
  expect_length(events,1L)
  expect_match(events[[1L]]$snapshot$render_line,"Refining bandwidth",fixed=TRUE)
  expect_match(events[[1L]]$snapshot$render_line,"attempting recovery",fixed=TRUE)
  expect_false(isTRUE(s$.np_progress_runtime$bandwidth_forward_active))
  expect_null(s$.np_progress_runtime$bandwidth_forward_state)
  s$.np_progress_abort(after,detail="original error")
  expect_identical(tail(events,1L)[[1L]]$event,"abort")
  expect_null(s$.np_progress_registry$active_id)
  events <- list()
  # The ended state must not print or reclaim ownership.
  s$.np_progress_bandwidth_nn_recovery()
  expect_length(events,0L)
})

test_that("recovery obeys ordinary message suppression and noninteractive policy", {
  for(kind in c("option-off","noninteractive","muffled")) {
    s <- nn_notice_scope(); events <- list()
    old <- options(np.messages=kind!="option-off")
    tryCatch({
      s$.np_progress_is_interactive <- function() kind!="noninteractive"
      s$.np_progress_resolve_message_muffling <-
        get(".np_progress_resolve_message_muffling",environment(.np_progress_begin))
      environment(s$.np_progress_resolve_message_muffling) <- s
      s$.np_progress_render_single_line <- function(snapshot,event)
        events[[length(events)+1L]] <<- snapshot
      s$.np_progress_runtime$bandwidth_state <- s$.np_progress_begin(
        "Bandwidth selection",domain="bandwidth")
      suppressMessages(s$.np_progress_bandwidth_nn_recovery())
      expect_length(events,0L)
      s$.np_progress_abort(s$.np_progress_runtime$bandwidth_state)
      expect_null(s$.np_progress_registry$active_id)
    },finally=options(old))
  }
})

test_that("recovery notices fit narrow consoles without claiming restoration", {
  for (prefix in c("[np]","[npRmpi]")) for (width in c(20L,40L,80L,180L)) {
    line <- .np_progress_bandwidth_notice_line(
      paste(prefix,"Bandwidth selection (multistart 10/10, iteration 640, elapsed 15.0s)"),
      prefix, c("NN bandwidth infeasible; attempting recovery","NN recovery"),
      width, fallback="recovering")
    expect_lte(nchar(line,type="width"),width)
    expect_match(line,"recover")
    expect_false(grepl("restored|finishing",line))
  }
})

test_that("R recovery notice neither changes probes nor labels skipped recovery", {
  s <- nn_notice_scope(); notices <- 0L; visited <- list()
  s$.np_progress_bandwidth_nn_recovery <- function() {notices <<- notices+1L; invisible(NULL)}
  result <- s$.np_nn_find_raw_valid_start(c(2,.3),1L,8L,function(p) {
    visited[[length(visited)+1L]] <<- p
    if(p[1L]<8) .Machine$double.xmax else 1.25
  })
  expect_identical(notices,1L)
  expect_identical(visited,list(c(4,.3),c(8,.3)))
  expect_identical(result,list(found=TRUE,point=c(8,.3),objective=1.25,evaluations=2L))
  empty <- s$.np_nn_find_raw_valid_start(c(12,.3),1L,8L,function(p) stop("not evaluated"))
  expect_false(empty$found)
  expect_identical(notices,1L)
  expect_error(s$.np_nn_find_raw_valid_start(c(2,.3),1L,8L,function(p) stop("original failure")),
               "original failure",fixed=TRUE)
  expect_identical(notices,2L)
})
