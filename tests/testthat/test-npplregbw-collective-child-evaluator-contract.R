run_npplreg_collective_child_evaluator_contract <- function() {
  close_mpi_slaves(force = TRUE)

  trace.file <- tempfile("npplreg-collective-trace-", fileext = ".log")
  old.trace.opt <- getOption("npRmpi.transport.trace.file")
  on.exit({
    options(npRmpi.transport.trace.file = old.trace.opt)
    unlink(trace.file)
  }, add = TRUE)

  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  ns <- asNamespace("npRmpi")
  options(npRmpi.transport.trace.file = trace.file)
  get(".npRmpi_bcast_cmd_expr", envir = ns, inherits = FALSE)(
    substitute(options(npRmpi.transport.trace.file = TF), list(TF = trace.file)),
    comm = 1L,
    caller.execute = TRUE
  )

  set.seed(20260408)
  n <- 24L
  xdat <- data.frame(x1 = runif(n), x2 = rnorm(n))
  zdat <- data.frame(z = sort(runif(n)))
  y <- 1 + xdat$x1 - 0.5 * xdat$x2 + sin(2 * pi * zdat$z) + rnorm(n, sd = 0.05)

  reg.args <- list(
    regtype = "lp",
    basis = "glp",
    degree = 1L,
    bernstein.basis = TRUE,
    bandwidth.compute = FALSE,
    bwmethod = "cv.ls",
    bwtype = "fixed"
  )
  outer.args <- reg.args[setdiff(names(reg.args), "bandwidth.compute")]
  bws <- matrix(0, nrow = 1L + ncol(xdat), ncol = ncol(zdat))
  degree.search <- list(
    engine = "nomad+powell",
    start.degree = 1L,
    candidates = list(c(0L, 1L)),
    lower = 0L,
    upper = 1L,
    basis = "glp",
    nobs = as.integer(n),
    start.user = FALSE,
    bernstein.basis = TRUE
  )

  state <- get(".npplregbw_nomad_prepare_state", envir = ns, inherits = FALSE)(
    xdat = xdat,
    ydat = y,
    zdat = zdat,
    bws = bws,
    reg.args = reg.args,
    outer.args = outer.args,
    degree.search = degree.search
  )

  bw.matrix <- matrix(
    c(0.7, 0.8, 0.9),
    nrow = length(state$child.templates),
    ncol = ncol(zdat)
  )

  serial <- get(".npplregbw_eval_child_payload_serial", envir = ns, inherits = FALSE)(
    zdat = zdat,
    reg.args = reg.args,
    degree.search = degree.search,
    child.responses = state$child.responses,
    child.templates = state$child.templates,
    child.setup = state$child.setup,
    bw.matrix = bw.matrix,
    degree = 1L,
    penalty.multiplier = 10
  )

  collective.expr <- substitute(
    get(".npplregbw_eval_child_payload_collective", envir = asNamespace("npRmpi"), inherits = FALSE)(
      ZDAT,
      REGARGS,
      DEGREESEARCH,
      CHILDRESPONSES,
      CHILDTEMPLATES,
      CHILDSETUP,
      BWMATRIX,
      DEGREE,
      PENALTY,
      COMM,
      EVALID
    ),
    list(
      ZDAT = zdat,
      REGARGS = reg.args,
      DEGREESEARCH = degree.search,
      CHILDRESPONSES = state$child.responses,
      CHILDTEMPLATES = state$child.templates,
      CHILDSETUP = state$child.setup,
      BWMATRIX = bw.matrix,
      DEGREE = 1L,
      PENALTY = 10,
      COMM = 1L,
      EVALID = 1L
    )
  )

  collective <- get(".npRmpi_bcast_cmd_expr", envir = ns, inherits = FALSE)(
    collective.expr,
    comm = 1L,
    caller.execute = TRUE
  )

  expect_equal(collective$objective, serial$objective, tolerance = 1e-10)
  expect_equal(collective$num.feval, serial$num.feval, tolerance = 1e-10)
  expect_equal(collective$num.feval.fast, serial$num.feval.fast, tolerance = 1e-10)

  trace.lines <- readLines(trace.file, warn = FALSE)
  start.lines <- grep(
    "role=npplreg.nomad\tevent=collective.eval.start",
    trace.lines,
    value = TRUE
  )
  done.lines <- grep(
    "role=npplreg.nomad\tevent=collective.eval.done",
    trace.lines,
    value = TRUE
  )

  expect_true(length(start.lines) >= 2L)
  expect_true(length(done.lines) >= 2L)
  expect_true(any(grepl("rank=0", start.lines, fixed = TRUE)))
  expect_true(any(grepl("rank=1", start.lines, fixed = TRUE)))
  expect_true(any(grepl("child_indices=1($|\t)", start.lines)))
  expect_true(any(grepl("child_indices=2,3($|\t)", start.lines)))
}

test_that("npplreg collective child evaluator matches the serial oracle at a fixed point", {
  skip_on_cran()

  if (.mpi_suite_pool_owned()) {
    marker <- "NPPLREGBW_COLLECTIVE_CHILD_ISOLATED_OK"
    test.path <- normalizePath(
      testthat::test_path("test-npplregbw-collective-child-evaluator-contract.R"),
      mustWork = TRUE
    )
    res <- npRmpi_run_isolated_contract(
      lines = c(
        "suppressPackageStartupMessages(library(testthat))",
        "suppressPackageStartupMessages(library(npRmpi))",
        sprintf(
          paste0(
            "testthat::test_file(%s, package = 'npRmpi', ",
            "load_package = 'none', stop_on_failure = TRUE)"
          ),
          dQuote(test.path)
        ),
        sprintf("cat('%s\\\\n')", marker)
      ),
      timeout = 180L,
      marker = marker
    )

    expect_false(is.null(res), info = "exact npRmpi install unavailable to isolated test")
    expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
    expect_true(res$witnessed, info = paste(res$output, collapse = "\n"))
  } else {
    run_npplreg_collective_child_evaluator_contract()
  }
})
