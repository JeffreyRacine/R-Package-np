test_that("SPMD local execution retains complete arbitrary application errors", {
  registry <- getFromNamespace(".npRmpi_spmd_registry", "npRmpi")
  opcode <- "unit.original.condition"
  old.seq <- getOption("npRmpi.spmd.seq_id", 0L)
  on.exit({options(npRmpi.spmd.seq_id = old.seq)
    if (exists(opcode, envir = registry, inherits = FALSE))
      rm(list = opcode, envir = registry)
  }, add = TRUE)
  original <- structure(simpleError("typed application failure", quote(f(x))),
    class = c("unit_typed_error", "error", "condition"),
    detail = list(code = 17L, field = "y"))
  assign(opcode, function(payload, envelope) stop(payload), envir = registry)
  envelope <- getFromNamespace(".npRmpi_spmd_make_envelope", "npRmpi")(
    opcode = opcode, timeout_class = "unit")
  captured <- getFromNamespace(".npRmpi_spmd_try_execute_local", "npRmpi")(
    envelope, original)
  expect_false(captured$ok)
  expect_identical(captured$condition, original)
  expect_identical(captured$error, conditionMessage(original))
  expect_identical(captured$ack$opcode, opcode)
  expect_identical(captured$ack$seq_id, envelope$seq_id)
})


test_that("command and argument evaluation preserve captured error objects", {
  original <- structure(simpleError("original", quote(f(x))),
    class = c("unit_typed_error", "error", "condition"), details = list(a = 1:3))
  envir <- new.env(parent = baseenv())
  envir$original <- original
  envir$f <- function(x) stop(original)
  envir$x <- 17L
  command <- getFromNamespace(".npRmpi_eval_scmd", "npRmpi")
  without <- getFromNamespace(".npRmpi_eval_without_dispatch", "npRmpi")
  resolve <- getFromNamespace(".npRmpi_autodispatch_resolve_owned_arg", "npRmpi")
  expect_identical(tryCatch(command(quote(f(x)), envir = envir), error = identity), original)
  expect_identical(tryCatch(command(envir$f, list(17L), envir = envir), error = identity), original)
  expect_identical(tryCatch(without(quote(f(x)), envir), error = identity), original)
  expect_identical(tryCatch(resolve(quote(f(x)), "owned", envir, character()), error = identity), original)
  # A returned error object is a value, not a signalled error.
  expect_identical(command(quote(original), envir = envir), original)
  expect_identical(command(quote(list(x, x + 1L)), envir = envir), list(17L,18L))
})
