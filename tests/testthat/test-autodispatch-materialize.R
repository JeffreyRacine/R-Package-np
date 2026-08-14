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

test_that("autodispatch remote references are reused only while value-current", {
  reset <- getFromNamespace(".npRmpi_lease_reset_local", "npRmpi")
  plan <- getFromNamespace(".npRmpi_lease_publication_plan", "npRmpi")
  prepare <- getFromNamespace(".npRmpi_lease_prepare_local", "npRmpi")
  commit <- getFromNamespace(".npRmpi_lease_commit_local", "npRmpi")
  tag <- getFromNamespace(".npRmpi_autodispatch_tag_result", "npRmpi")
  current <- getFromNamespace(".npRmpi_autodispatch_ref_is_current", "npRmpi")
  reset()
  on.exit(reset(), add = TRUE)

  value <- structure(list(bws = c(1, 2), degree = c(1L, 1L)), class = "rbandwidth")
  publication <- plan(quote(npregbw(xdat = x, ydat = y)))
  prepare(value, publication)
  commit(publication)
  value <- tag(value, publication = publication)
  expect_true(current(value))
  value$bws[[1L]] <- 3
  expect_false(current(value))
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

test_that("autodispatch materializes forwarded formals and named dots from one owner frame", {
  materialize <- getFromNamespace(".npRmpi_autodispatch_materialize_call", "npRmpi")

  owner <- function(dat, bws, ...) {
    mc <- match.call()
    mc[[1L]] <- as.name("npudistbw")
    materialize(mc = mc, caller_env = environment(), comm = 1L)
  }
  forward <- function(...) owner(...)

  dat.value <- c(2, 7, 1, 8)
  bw.value <- 0.35
  nmulti <- 99L
  nomad.nmulti <- 77L
  prepared <- forward(
    dat = dat.value,
    bws = bw.value,
    nmulti = 1L,
    nomad.nmulti = 0L,
    proper.control = NULL
  )

  expect_identical(prepared$tmpvals[[as.character(prepared$call$dat)]], dat.value)
  expect_identical(prepared$tmpvals[[as.character(prepared$call$bws)]], bw.value)
  expect_identical(prepared$tmpvals[[as.character(prepared$call$nmulti)]], 1L)
  expect_identical(prepared$tmpvals[[as.character(prepared$call$nomad.nmulti)]], 0L)
  expect_true("proper.control" %in% names(as.list(prepared$call)))
  expect_null(prepared$call$proper.control)
})

test_that("autodispatch does not discover unresolved arguments in dynamic frames", {
  materialize <- getFromNamespace(".npRmpi_autodispatch_materialize_call", "npRmpi")

  leaked.control <- 4L
  owner <- new.env(parent = baseenv())
  expect_error(
    materialize(
      mc = quote(npudistbw(nmulti = leaked.control)),
      caller_env = owner,
      comm = 1L
    ),
    "leaked.control"
  )
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

test_that("autodispatch materialization resolves estimator uncertainty controls", {
  materialize <- getFromNamespace(".npRmpi_autodispatch_materialize_call", "npRmpi")

  bws <- structure(list(tag = "dummy"), class = "lsqregressionbandwidth")
  gradients <- FALSE
  residuals <- TRUE
  se <- FALSE

  mc <- quote(nplsqreg(
    bws = bws,
    gradients = gradients,
    residuals = residuals,
    se = se
  ))
  prepared <- materialize(mc = mc, caller_env = environment(), comm = 1L)

  gradients.ref <- as.character(prepared$call$gradients)
  residuals.ref <- as.character(prepared$call$residuals)
  se.ref <- as.character(prepared$call$se)

  expect_identical(prepared$tmpvals[[gradients.ref]], FALSE)
  expect_identical(prepared$tmpvals[[residuals.ref]], TRUE)
  expect_identical(prepared$tmpvals[[se.ref]], FALSE)
  expect_false(identical(prepared$call$se, as.name("se")))
})

test_that("registered autodispatch call formals have an explicit transport owner", {
  ns <- asNamespace("npRmpi")
  targets <- getFromNamespace(".npRmpi_autodispatch_target_args", "npRmpi")()

  character_constants <- function(expr) {
    if (is.character(expr))
      return(expr)
    if (!is.recursive(expr))
      return(character(0))
    unlist(lapply(as.list(expr), character_constants), use.names = FALSE)
  }

  helper.names <- ls(ns, all.names = TRUE)
  helper.names <- helper.names[
    grepl("^\\.npRmpi_autodispatch_is_.*_core$", helper.names)
  ]
  call.heads <- unique(unlist(lapply(helper.names, function(nm) {
    character_constants(body(get(nm, envir = ns, inherits = FALSE)))
  }), use.names = FALSE))
  call.heads <- unique(c(
    call.heads,
    "npregbw", "npscoefbw", "npplregbw", "npindexbw",
    "npudensbw", "npudistbw", "npcdensbw", "npcdistbw"
  ))
  call.heads <- call.heads[vapply(call.heads, exists, logical(1),
                                  envir = ns, mode = "function",
                                  inherits = FALSE)]

  formal.pairs <- do.call(rbind, lapply(sort(call.heads), function(head) {
    args <- setdiff(names(formals(get(head, envir = ns, inherits = FALSE))), "...")
    if (!length(args))
      return(NULL)
    data.frame(head = rep(head, length(args)), arg = args,
               stringsAsFactors = FALSE)
  }))

  # `subset` must remain language for model-frame evaluation.  The internal
  # formula-method `call` formal is deliberately reconstructed by its owner.
  language.or.derived <- c("subset", "call")
  unowned <- formal.pairs[
    !formal.pairs$arg %in% c(targets, language.or.derived), , drop = FALSE
  ]

  expect_equal(unowned, formal.pairs[FALSE, , drop = FALSE],
               info = paste(capture.output(print(unowned)), collapse = "\n"))
})

test_that("ordinary autodispatch controls are evaluated in the caller", {
  materialize <- getFromNamespace(".npRmpi_autodispatch_materialize_call", "npRmpi")

  probabilities <- TRUE
  level <- "high"
  pivot <- FALSE
  bw.x <- 0.25
  maxiter <- 17L
  mc <- quote(npconmode(
    probabilities = probabilities,
    level = level,
    pivot = pivot,
    bw.x = bw.x,
    maxiter = maxiter
  ))
  prepared <- materialize(mc = mc, caller_env = environment(), comm = 1L)

  for (nm in c("probabilities", "level", "pivot", "bw.x", "maxiter")) {
    ref <- as.character(prepared$call[[nm]])
    expect_identical(prepared$tmpvals[[ref]], get(nm, inherits = FALSE))
  }
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
    powell.remin = FALSE,
    itmax = 1L
  ))

  prepared <- materialize(mc = mc, caller_env = environment(), comm = 1L)
  call.args <- as.list(prepared$call)

  expect_true("gydat" %in% names(call.args))
  expect_true(is.null(call.args$gydat))

  bws.ref <- as.character(prepared$call$bws)
  compute.ref <- as.character(prepared$call$bandwidth.compute)
  nmulti.ref <- as.character(prepared$call$nmulti)
  remin.ref <- as.character(prepared$call$powell.remin)
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
  expect_identical(prepared$tmpvals[[as.character(prepared$call$bwmethod)]], "cv.ls")
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
