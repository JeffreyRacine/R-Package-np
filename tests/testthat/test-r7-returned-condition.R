test_that("frame evaluation distinguishes returned values from raised conditions", {
  env <- new.env(parent = baseenv())
  env$err <- structure(list(message = "r7 probe", call = NULL, marker = 61L),
                       class = c("r7_probe", "error", "condition"))
  env$attempts <- 0L
  for (value in list(env$err, structure("returned text", class = "try-error"),
                    NULL, list(ok = FALSE, error = env$err), 42L)) {
    env$value <- value
    for (expr in list(quote(identity(value)), quote(value))) {
      out <- .np_try_eval_in_frames(expr, eval_env = env, search_frames = FALSE)
      expect_true(out$ok)
      expect_identical(out$value, value)
      expect_null(out$error)
    }
    out <- .np_try_eval_in_frames(quote(identity(value)),
      eval_env = list(value = value), enclos = baseenv(), search_frames = FALSE)
    expect_true(out$ok)
    expect_identical(out$value, value)
  }
  caught <- .np_try_eval_in_frames(
    quote(tryCatch(stop(err), error = identity)), env, search_frames = FALSE)
  expect_true(caught$ok)
  expect_identical(caught$value, env$err)
  out <- .np_try_eval_in_frames(
    quote({ attempts <- attempts + 1L; stop(err) }), env, search_frames = FALSE)
  expect_false(out$ok)
  expect_identical(out$error, env$err)
  expect_identical(env$attempts, 1L)
  env$x <- 9L
  x <- 99L
  expect_identical(.np_try_eval_in_frames(quote(x + 1L), env)$value, 10L)
})

test_that("worker command values do not impersonate caught exceptions", {
  loop <- .npRmpi_worker_loop
  env <- new.env(parent = environment(loop))
  environment(loop) <- env
  for (value in list(simpleError("returned"),
                    structure("returned try", class = "try-error"), NULL)) {
    count <- 0L
    env$mpi.bcast.cmd <- function(...) {
      count <<- count + 1L
      if (count == 1L) quote(identity(1L)) else "kaerb"
    }
    env$.npRmpi_eval_scmd <- function(...) value
    expect_silent(loop())
    expect_identical(count, 2L)
  }
  env$mpi.bcast.cmd <- function(...) quote(stop("raised"))
  env$.npRmpi_eval_scmd <- function(...) stop("raised", call. = FALSE)
  env$mpi.comm.rank <- function(...) 1L
  expect_error(loop(), "raised")
})
