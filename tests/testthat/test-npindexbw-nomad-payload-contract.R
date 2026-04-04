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
        engine = "nomad+powell"
      ),
      nomad.inner.nmulti = 0L
    )
  )

  expect_equal(result$payload$fval, 10)
  expect_equal(result$payload$num.feval, 8)
  expect_equal(result$payload$num.feval.fast, 8)
  expect_equal(result$objective, 10)
})
