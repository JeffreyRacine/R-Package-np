test_that("raw-state replay refuses a remote custom generator before advance", {
  en <- new.env(parent=environment(.npRmpi_bootstrap_attach_rng_stream))
  f <- .npRmpi_bootstrap_attach_rng_stream
  environment(f) <- en
  en$RNGkind <- function() c("user-supplied", "Inversion", "Rejection")
  en$.npRmpi_has_active_slave_pool <- function() TRUE
  en$.npRmpi_master_only_mode <- function() FALSE
  set.seed(17); before <- .Random.seed; advanced <- FALSE
  tasks <- list(list(start=1L, bsz=3L))
  expect_error(f(tasks, function(task) { advanced <<- TRUE; runif(3) }),
    "cannot reproduce a user-supplied uniform RNG")
  expect_false(advanced)
  expect_identical(.Random.seed, before)
  # Built-in effective streams and explicit local routes retain their owner.
  en$RNGkind <- function() c("Mersenne-Twister", "Inversion", "Rejection")
  z <- f(tasks, function(task) runif(task$bsz))
  expect_identical(z[[1L]]$rng_state, before)
  expect_identical(attr(z, "rng_final_state"), .Random.seed)
  en$RNGkind <- function() c("user-supplied", "Inversion", "Rejection")
  en$.npRmpi_has_active_slave_pool <- function() FALSE
  en$.npRmpi_master_only_mode <- function() TRUE
  expect_no_error(f(tasks, function(task) runif(task$bsz)))
})


test_that("custom normal admission is limited to actual normal tasks", {
  f <- .npRmpi_bootstrap_attach_rng_stream
  en <- new.env(parent = environment(f)); environment(f) <- en
  en$RNGkind <- function() c("Mersenne-Twister", "user-supplied", "Rejection")
  queried <- 0L; advanced <- 0L
  en$.npRmpi_bootstrap_require_normal_rng <- function() {
    queried <<- queried + 1L
    stop("normal generator absent")
  }
  tasks <- list(list(start = 1L, bsz = 3L))
  advance <- function(task) { advanced <<- advanced + 1L; runif(task$bsz) }
  set.seed(71); before <- .Random.seed
  expect_error(f(tasks, advance, normal.draws = TRUE), "normal generator absent")
  expect_identical(queried, 1L); expect_identical(advanced, 0L)
  expect_identical(.Random.seed, before)
  uniform <- f(tasks, advance)
  expect_null(uniform[[1L]][["normal.kind", exact = TRUE]])
  expect_identical(queried, 1L)
  en$.npRmpi_bootstrap_require_normal_rng <- function() invisible(NULL)
  set.seed(71); normal <- f(tasks, advance, normal.draws = TRUE)
  expect_identical(normal[[1L]]$normal.kind, "user-supplied")
  expect_identical(normal[[1L]]$rng_state, before)
  normal[[1L]]$normal.kind <- NULL
  expect_identical(normal, uniform)
})

test_that("custom normal initialization precedes replay state assignment", {
  f <- .npRmpi_bootstrap_use_task_rng
  en <- new.env(parent = environment(f)); environment(f) <- en
  events <- character()
  en$RNGkind <- function(normal.kind) events <<- c(events, normal.kind)
  en$assign <- function(x, value, envir) events <<- c(events, "seed")
  f(list(normal.kind = "user-supplied", rng_state = 1:3))
  expect_identical(events, c("user-supplied", "seed"))
  events <- character(); f(list(rng_state = 1:3))
  expect_identical(events, "seed")
})
