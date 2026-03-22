test_that("autodispatch materialization preserves explicit argument expressions", {
  materialize <- getFromNamespace(".npRmpi_autodispatch_materialize_call", "npRmpi")

  bws <- list(bw = c(0.11, 0.22))
  bw.j <- bws$bw[1L]
  x.marginal <- c(1, 2, 3)

  mc <- quote(npudistbw(dat = x.marginal, bws = bw.j, bandwidth.compute = FALSE))
  prepared <- materialize(mc = mc, caller_env = environment(), comm = 1L)

  dat.ref <- as.character(prepared$call$dat)
  bws.ref <- as.character(prepared$call$bws)

  expect_identical(prepared$tmpvals[[dat.ref]], x.marginal)
  expect_identical(prepared$tmpvals[[bws.ref]], bw.j)
  expect_false(identical(prepared$tmpvals[[bws.ref]], bws))
})

test_that("autodispatch materialization resolves ..n placeholders by argument name", {
  materialize <- getFromNamespace(".npRmpi_autodispatch_materialize_call", "npRmpi")

  dat <- c(3, 1, 4)
  bw.j <- 0.5
  mc <- as.call(list(
    as.name("npudistbw"),
    dat = as.name("..1"),
    bws = quote(bw.j),
    bandwidth.compute = FALSE
  ))

  prepared <- materialize(mc = mc, caller_env = environment(), comm = 1L)
  dat.ref <- as.character(prepared$call$dat)
  bws.ref <- as.character(prepared$call$bws)

  expect_identical(prepared$tmpvals[[dat.ref]], dat)
  expect_identical(prepared$tmpvals[[bws.ref]], bw.j)
})

test_that("autodispatch materialization evaluates proper arguments eagerly", {
  materialize <- getFromNamespace(".npRmpi_autodispatch_materialize_call", "npRmpi")

  bws <- structure(list(x = 1), class = "conbandwidth")
  proper.flag <- TRUE
  proper.method <- "project"
  proper.control <- list(mode = "slice", slice.grid.size = 21L)

  mc <- quote(npcdens(
    bws = bws,
    proper = proper.flag,
    proper.method = proper.method,
    proper.control = proper.control
  ))
  prepared <- materialize(mc = mc, caller_env = environment(), comm = 1L)

  proper.ref <- as.character(prepared$call$proper)
  method.ref <- as.character(prepared$call$proper.method)
  control.ref <- as.character(prepared$call$proper.control)

  expect_identical(prepared$tmpvals[[proper.ref]], proper.flag)
  expect_identical(prepared$tmpvals[[method.ref]], proper.method)
  expect_identical(prepared$tmpvals[[control.ref]], proper.control)
})
