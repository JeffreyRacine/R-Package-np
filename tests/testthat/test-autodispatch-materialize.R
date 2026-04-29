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

test_that("autodispatch materialization preserves explicit NULL arguments without shifting later args", {
  materialize <- getFromNamespace(".npRmpi_autodispatch_materialize_call", "npRmpi")

  xdat <- data.frame(x = 1:3)
  ydat <- c(1, 2, 3)
  gydat <- NULL
  bws <- structure(list(tag = "dummy"), class = "condbandwidth")

  mc <- quote(npcdistbw(
    xdat = xdat,
    ydat = ydat,
    gydat = gydat,
    bws = bws,
    bandwidth.compute = TRUE,
    nmulti = 1L,
    remin = FALSE,
    itmax = 1L
  ))

  prepared <- materialize(mc = mc, caller_env = environment(), comm = 1L)
  call.args <- as.list(prepared$call)

  expect_true("gydat" %in% names(call.args))
  expect_true(is.null(call.args$gydat))

  bws.ref <- as.character(prepared$call$bws)
  compute.ref <- as.character(prepared$call$bandwidth.compute)
  nmulti.ref <- as.character(prepared$call$nmulti)
  remin.ref <- as.character(prepared$call$remin)
  itmax.ref <- as.character(prepared$call$itmax)

  expect_false(identical(bws.ref, "bws"))
  expect_identical(prepared$tmpvals[[compute.ref]], TRUE)
  expect_identical(prepared$tmpvals[[nmulti.ref]], 1L)
  expect_identical(prepared$tmpvals[[remin.ref]], FALSE)
  expect_identical(prepared$tmpvals[[itmax.ref]], 1L)
})

test_that("autodispatch helpers leave calls without deferred dots unchanged", {
  genericize <- getFromNamespace(".npRmpi_autodispatch_as_generic_call", "npRmpi")
  materialize <- getFromNamespace(".npRmpi_autodispatch_materialize_call", "npRmpi")

  xdat <- data.frame(x = 1:3)
  ydat <- c(1, 2, 3)
  bws <- 0
  mc <- quote(npregbw.NULL(
    xdat = xdat,
    ydat = ydat,
    bws = bws,
    regtype = "lc",
    bwmethod = "cv.ls"
  ))

  rebuilt <- genericize("npregbw", mc)
  expect_identical(rebuilt, quote(npregbw(
    xdat = xdat,
    ydat = ydat,
    bws = bws,
    regtype = "lc",
    bwmethod = "cv.ls"
  )))

  prepared <- materialize(mc = rebuilt, caller_env = environment(), comm = 1L)
  call.args <- as.list(prepared$call)[-1L]

  expect_false("..." %in% names(call.args))
  expect_true(all(c("xdat", "ydat", "bws", "regtype", "bwmethod") %in% names(call.args)))
  expect_identical(prepared$tmpvals[[as.character(prepared$call$regtype)]], "lc")
  expect_identical(prepared$call$bwmethod, "cv.ls")
})

test_that("autodispatch generic-call rebuild expands deferred dots from method match.call", {
  genericize <- getFromNamespace(".npRmpi_autodispatch_as_generic_call", "npRmpi")

  mc <- quote(npregbw.NULL(
    xdat = xdat,
    ydat = ydat,
    bws = bws,
    ... = pairlist(regtype = "ll", bwmethod = "cv.aic")
  ))

  rebuilt <- genericize("npregbw", mc)
  rebuilt.args <- as.list(rebuilt)[-1L]

  expect_identical(as.character(rebuilt[[1L]]), "npregbw")
  expect_true(all(c("xdat", "ydat", "bws", "regtype", "bwmethod") %in% names(rebuilt.args)))
  expect_false("..." %in% names(rebuilt.args))
  expect_identical(rebuilt.args$regtype, "ll")
  expect_identical(rebuilt.args$bwmethod, "cv.aic")
})

test_that("autodispatch materialization expands deferred dots before shipping attach calls", {
  materialize <- getFromNamespace(".npRmpi_autodispatch_materialize_call", "npRmpi")

  xdat <- data.frame(x1 = 1:3, x2 = 4:6)
  ydat <- c(0, 1, 0)
  bws <- 0
  mc <- quote(npindexbw.NULL(
    xdat = xdat,
    ydat = ydat,
    bws = bws,
    ... = pairlist(method = "kleinspady", regtype = "ll")
  ))

  prepared <- materialize(mc = mc, caller_env = environment(), comm = 1L)
  call.args <- as.list(prepared$call)[-1L]

  expect_true(all(c("xdat", "ydat", "bws", "method", "regtype") %in% names(call.args)))
  expect_false("..." %in% names(call.args))

  method.ref <- as.character(prepared$call$method)
  regtype.ref <- as.character(prepared$call$regtype)

  expect_identical(prepared$tmpvals[[method.ref]], "kleinspady")
  expect_identical(prepared$tmpvals[[regtype.ref]], "ll")
})
