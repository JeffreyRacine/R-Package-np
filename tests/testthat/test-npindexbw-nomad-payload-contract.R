with_nprmpi_npindex_nomad_payload_bindings <- function(bindings, code) {
  code <- substitute(code)
  ns <- asNamespace("npRmpi")
  old <- lapply(names(bindings), function(name) get(name, envir = ns, inherits = FALSE))
  names(old) <- names(bindings)

  for (name in names(bindings)) {
    was_locked <- bindingIsLocked(name, ns)
    if (was_locked)
      unlockBinding(name, ns)
    assign(name, bindings[[name]], envir = ns)
    if (was_locked)
      lockBinding(name, ns)
  }

  on.exit({
    for (name in names(old)) {
      was_locked <- bindingIsLocked(name, ns)
      if (was_locked)
        unlockBinding(name, ns)
      assign(name, old[[name]], envir = ns)
      if (was_locked)
        lockBinding(name, ns)
    }
  }, add = TRUE)

  eval(code, envir = parent.frame())
}

test_that("npindexbw retains Powell evaluation counts when Powell does not improve", {
  ns <- asNamespace("npRmpi")

  result <- with_nprmpi_npindex_nomad_payload_bindings(
    list(
      .npindexbw_build_sibandwidth = function(xdat, ydat, bws, template, bandwidth.compute, reg.args) {
        list(
          beta = c(1, 0.5),
          bw = 0.3,
          type = "fixed",
          method = "ichimura",
          ckertype = "gaussian",
          ckerorder = 2L,
          ckerbound = "none",
          ckerlb = NA_real_,
          ckerub = NA_real_,
          ynames = "y"
        )
      },
      .npindexbw_eval_objective = function(param, xmat, ydat, bws, spec) {
        list(objective = 10, num.feval.fast = 3L)
      },
      .np_nomad_baseline_note = function(degree) invisible(NULL),
      .np_nomad_search = function(engine,
                                  baseline_record,
                                  start_degree,
                                  x0,
                                  bbin,
                                  lb,
                                  ub,
                                  eval_fun,
                                  build_payload,
                                  direction,
                                  objective_name,
                                  nmulti,
                                  nomad.inner.nmulti,
                                  random.seed,
                                  degree_spec,
                                  progress_label) {
        eval_fun(x0)
        build_payload(
          x0,
          best_record = list(degree = start_degree, objective = 10, num.feval = 1L),
          solution = list(bbe = 1L),
          interrupted = FALSE
        )
      },
      .np_nomad_with_powell_progress = function(degree, expr) expr,
      .npindexbw_run_fixed_degree = function(xdat, ydat, bws, template, reg.args, opt.args) {
        list(
          method = "ichimura",
          fval = 11,
          num.feval = 7L,
          num.feval.fast = 5L
        )
      },
      npindexbw.sibandwidth = function(xdat, ydat, bws, bandwidth.compute, ...) {
        bws
      }
    ),
    get(".npindexbw_nomad_search", envir = ns, inherits = FALSE)(
      xdat = data.frame(x1 = c(0.1, 0.2), x2 = c(0.3, 0.4)),
      ydat = c(1, 2),
      bws = c(1, 0.5, 0.3),
      template = list(method = "ichimura"),
      reg.args = list(
        method = "ichimura",
        regtype = "lp",
        basis.engine = "tensor"
      ),
      opt.args = list(nmulti = 1L),
      degree.search = list(
        verify = FALSE,
        start.degree = 1L,
        bernstein.basis = TRUE,
        lower = 0L,
        upper = 2L,
        candidates = list(0:2),
        engine = "nomad+powell",
        basis = "glp",
        nobs = 2L,
        start.user = FALSE
      ),
      nomad.inner.nmulti = 0L
    )
  )

  expect_equal(result$payload$fval, 10)
  expect_equal(result$payload$num.feval, 8)
  expect_equal(result$payload$num.feval.fast, 8)
  expect_equal(result$objective, 10)
})

test_that("npindexbw fixed NOMAD route normalizes internal h starts but keeps public restart starts raw", {
  ns <- asNamespace("npRmpi")
  xdat <- data.frame(x1 = c(-0.3, 0.1, 0.4, 0.8), x2 = c(0.5, -0.1, 0.2, -0.4))
  ydat <- c(0.2, -0.1, 0.8, 1.1)
  baseline.bws <- list(
    beta = c(1, 0),
    bw = 0,
    type = "fixed",
    method = "ichimura",
    ckertype = "gaussian",
    ckerorder = 2L,
    ckerbound = "none",
    ckerlb = NA_real_,
    ckerub = NA_real_,
    ynames = "y"
  )
  degree.search <- list(
    verify = FALSE,
    start.degree = 1L,
    bernstein.basis = TRUE,
    lower = 0L,
    upper = 2L,
    candidates = list(0:2),
    engine = "nomad+powell",
    basis = "glp",
    nobs = nrow(xdat),
    start.user = FALSE
  )
  setup <- get(".npindexbw_nomad_fixed_start_setup", envir = ns, inherits = FALSE)(
    xmat = as.matrix(xdat),
    ydat = ydat,
    baseline.bws = baseline.bws,
    degree.search = degree.search,
    nmulti = 1L,
    random.seed = 42L
  )
  captured <- NULL

  result <- with_nprmpi_npindex_nomad_payload_bindings(
    list(
      .npindexbw_build_sibandwidth = function(xdat, ydat, bws, template, bandwidth.compute, reg.args) {
        baseline.bws
      },
      .np_nomad_baseline_note = function(degree) invisible(NULL),
      .np_nomad_search = function(engine,
                                  baseline_record,
                                  start_degree,
                                  x0,
                                  bbin,
                                  lb,
                                  ub,
                                  eval_fun,
                                  build_payload,
                                  direction,
                                  objective_name,
                                  nmulti,
                                  nomad.inner.nmulti,
                                  random.seed,
                                  degree_spec,
                                  ...) {
        captured <<- list(
          x0 = as.numeric(x0),
          lb = as.numeric(lb),
          ub = as.numeric(ub)
        )
        list(
          method = engine,
          direction = direction,
          verify = FALSE,
          completed = TRUE,
          certified = FALSE,
          interrupted = FALSE,
          baseline = list(degree = start_degree, objective = 10),
          best = list(degree = start_degree, objective = 10, num.feval = 1L),
          best_payload = list(fval = 10, num.feval = 1L, num.feval.fast = 0L, method = "ichimura"),
          best_point = as.numeric(x0),
          n.unique = 1L,
          n.visits = 1L,
          n.cached = 0L,
          nomad.time = 1,
          powell.time = 0,
          optim.time = 1,
          grid.size = NA_integer_,
          best.restart = 1L,
          restart.starts = list(as.numeric(x0)),
          restart.degree.starts = list(as.integer(x0[length(x0)])),
          restart.bandwidth.starts = list(as.numeric(x0[seq_len(length(x0) - 1L)])),
          restart.start.info = list(user_supplied_start = FALSE),
          restart.results = list(list(
            restart = 1L,
            start = as.numeric(x0),
            degree.start = start_degree,
            objective = 10,
            bbe = 1L,
            iterations = 0L,
            solution = as.numeric(x0)
          )),
          trace = data.frame()
        )
      }
    ),
    get(".npindexbw_nomad_search", envir = ns, inherits = FALSE)(
      xdat = xdat,
      ydat = ydat,
      bws = c(0, 0, 0),
      template = list(method = "ichimura"),
      reg.args = list(
        method = "ichimura",
        regtype = "lp",
        basis.engine = "tensor"
      ),
      opt.args = list(nmulti = 1L, random.seed = 42L),
      degree.search = degree.search,
      nomad.inner.nmulti = 0L
    )
  )

  expect_equal(captured$x0, as.numeric(setup$start_matrix.point[1L, ]), tolerance = 1e-12)
  expect_equal(result$restart.starts[[1]], as.numeric(setup$start_matrix.raw[1L, ]), tolerance = 1e-12)
  expect_equal(result$restart.bandwidth.starts[[1]], as.numeric(setup$start_matrix.raw[1L, seq_len(ncol(setup$start_matrix.raw) - 1L)]), tolerance = 1e-12)
  expect_equal(result$restart.results[[1]]$start, as.numeric(setup$start_matrix.raw[1L, ]), tolerance = 1e-12)
  expect_equal(result$restart.results[[1]]$solution, as.numeric(setup$start_matrix.raw[1L, ]), tolerance = 1e-12)
})

test_that("npindexbw fixed NOMAD route preserves explicit user fixed starts", {
  ns <- asNamespace("npRmpi")
  xdat <- data.frame(x1 = c(-0.3, 0.1, 0.4, 0.8), x2 = c(0.5, -0.1, 0.2, -0.4))
  ydat <- c(0.2, -0.1, 0.8, 1.1)
  baseline.bws <- list(
    beta = c(1, 1.75),
    bw = 0.25,
    type = "fixed",
    method = "ichimura",
    ckertype = "gaussian",
    ckerorder = 2L,
    ckerbound = "none",
    ckerlb = NA_real_,
    ckerub = NA_real_,
    ynames = "y"
  )
  degree.search <- list(
    verify = FALSE,
    start.degree = 1L,
    bernstein.basis = TRUE,
    lower = 0L,
    upper = 2L,
    candidates = list(0:2),
    engine = "nomad+powell",
    basis = "glp",
    nobs = nrow(xdat),
    start.user = FALSE
  )
  setup <- get(".npindexbw_nomad_fixed_start_setup", envir = ns, inherits = FALSE)(
    xmat = as.matrix(xdat),
    ydat = ydat,
    baseline.bws = baseline.bws,
    degree.search = degree.search,
    nmulti = 1L,
    random.seed = 42L
  )
  captured <- NULL

  result <- with_nprmpi_npindex_nomad_payload_bindings(
    list(
      .npindexbw_build_sibandwidth = function(xdat, ydat, bws, template, bandwidth.compute, reg.args) {
        baseline.bws
      },
      .np_nomad_baseline_note = function(degree) invisible(NULL),
      .np_nomad_search = function(engine,
                                  baseline_record,
                                  start_degree,
                                  x0,
                                  bbin,
                                  lb,
                                  ub,
                                  eval_fun,
                                  build_payload,
                                  direction,
                                  objective_name,
                                  nmulti,
                                  nomad.inner.nmulti,
                                  random.seed,
                                  degree_spec,
                                  ...) {
        captured <<- list(
          x0 = as.numeric(x0),
          lb = as.numeric(lb),
          ub = as.numeric(ub)
        )
        list(
          method = engine,
          direction = direction,
          verify = FALSE,
          completed = TRUE,
          certified = FALSE,
          interrupted = FALSE,
          baseline = list(degree = start_degree, objective = 10),
          best = list(degree = start_degree, objective = 10, num.feval = 1L),
          best_payload = list(fval = 10, num.feval = 1L, num.feval.fast = 0L, method = "ichimura"),
          best_point = as.numeric(x0),
          n.unique = 1L,
          n.visits = 1L,
          n.cached = 0L,
          nomad.time = 1,
          powell.time = 0,
          optim.time = 1,
          grid.size = NA_integer_,
          best.restart = 1L,
          restart.starts = list(as.numeric(x0)),
          restart.degree.starts = list(as.integer(x0[length(x0)])),
          restart.bandwidth.starts = list(as.numeric(x0[seq_len(length(x0) - 1L)])),
          restart.start.info = list(user_supplied_start = FALSE),
          restart.results = list(list(
            restart = 1L,
            start = as.numeric(x0),
            degree.start = start_degree,
            objective = 10,
            bbe = 1L,
            iterations = 0L,
            solution = as.numeric(x0)
          )),
          trace = data.frame()
        )
      }
    ),
    get(".npindexbw_nomad_search", envir = ns, inherits = FALSE)(
      xdat = xdat,
      ydat = ydat,
      bws = c(0, 0, 0),
      template = list(method = "ichimura"),
      reg.args = list(
        method = "ichimura",
        regtype = "lp",
        basis.engine = "tensor"
      ),
      opt.args = list(nmulti = 1L, random.seed = 42L),
      degree.search = degree.search,
      nomad.inner.nmulti = 0L
    )
  )

  expect_equal(setup$start_matrix.raw[1L, 1L], 1.75, tolerance = 1e-12)
  expect_equal(setup$start_matrix.raw[1L, 2L], 0.25, tolerance = 1e-12)
  expect_equal(captured$x0, as.numeric(setup$start_matrix.point[1L, ]), tolerance = 1e-12)
  expect_equal(result$restart.starts[[1]], as.numeric(setup$start_matrix.raw[1L, ]), tolerance = 1e-12)
  expect_equal(result$restart.results[[1]]$start, as.numeric(setup$start_matrix.raw[1L, ]), tolerance = 1e-12)
})
