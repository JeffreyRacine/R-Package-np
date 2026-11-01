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
