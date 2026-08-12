npcdens_y_side_source <- function(file) {
  candidates <- c(
    test_path("..", "..", "src", file),
    test_path("..", "..", "..", "src", file),
    file.path(Sys.getenv("R_PACKAGE_DIR", ""), "src", file),
    file.path(Sys.getenv("R_PACKAGE_SOURCE", ""), "src", file),
    file.path(getwd(), "src", file),
    file.path(getwd(), "..", "src", file)
  )
  candidates <- unique(candidates[nzchar(candidates)])
  hits <- candidates[file.exists(candidates)]
  if (!length(hits))
    return(NULL)
  paste(readLines(hits[[1L]], warn = FALSE), collapse = "\n")
}

npcdens_y_side_fixture <- function() {
  set.seed(20260726)
  n <- 48L
  x <- data.frame(
    x1 = runif(n),
    x2 = runif(n),
    u = factor(rep(c("a", "b"), length.out = n))
  )
  y <- data.frame(
    y = factor(ifelse(
      x$x1 + x$x2 + rnorm(n, sd = 0.15) > 1,
      "yes",
      "no"
    ))
  )
  list(x = x, y = y)
}

npcdens_y_side_oracle <- function(x, y, bw) {
  xweights <- npksum(
    txdat = x,
    exdat = x,
    bws = bw$xbw,
    bwtype = bw$type,
    ckertype = bw$cxkertype,
    ckerorder = bw$cxkerorder,
    ukertype = bw$uxkertype,
    okertype = bw$oxkertype,
    return.kernel.weights = TRUE
  )$kw
  yweights <- npksum(
    txdat = y,
    exdat = y,
    bws = bw$ybw,
    bwtype = bw$type,
    ukertype = bw$uykertype,
    okertype = bw$oykertype,
    return.kernel.weights = TRUE
  )$kw
  diag(xweights) <- 0
  diag(yweights) <- 0
  sum(log(colSums(xweights * yweights) / colSums(xweights)))
}

test_that("conditional CVML Y-side ownership follows the selected stream", {
  source <- npcdens_y_side_source("np.c")
  skip_if(is.null(source), "source file src/np.c unavailable")

  expect_true(grepl(
    "static int np_conditional_density_objective_needs_y_side(",
    source,
    fixed = TRUE
  ))
  expect_true(grepl(
    "np_conditional_density_cvml_stream_engine_supported());",
    source,
    fixed = TRUE
  ))
  expect_equal(
    lengths(regmatches(
      source,
      gregexpr(
        "need_y_side = np_conditional_density_objective_needs_y_side(ibwmfunc);",
        source,
        fixed = TRUE
      )
    )),
    1L
  )
  expect_false(grepl(
    "need_y_side = (ibwmfunc == CBWM_CVLS) || ((ibwmfunc == CBWM_CVML) && (np_lp_engine_extern == NP_LP_ENGINE_GENERAL));",
    source,
    fixed = TRUE
  ))
})

test_that("categorical kernel rows reject missing category metadata", {
  source <- npcdens_y_side_source("jksum.c")
  skip_if(is.null(source), "source file src/jksum.c unavailable")

  start <- regexpr("kernel_weighted_sum_np_ctx_ex\\(", source)
  stop <- regexpr("int kernel_weighted_sum_np_ctx(", source,
                  fixed = TRUE)
  expect_gt(start, 0L)
  expect_gt(stop, start)
  owner <- substr(source, start, stop - 1L)
  expect_true(grepl(
    "((num_reg_unordered > 0) || (num_reg_ordered > 0))",
    owner,
    fixed = TRUE
  ))
  expect_true(grepl("(num_categories == NULL)", owner, fixed = TRUE))
})

test_that("native CVML degree buffers follow the resident degree-key contract", {
  source <- npcdens_y_side_source("np.c")
  skip_if(is.null(source), "source file src/np.c unavailable")

  start <- regexpr(
    "static int np_conditional_density_prepared_context_eval_native_raw(",
    source,
    fixed = TRUE
  )
  stop <- regexpr(
    "SEXP C_np_density_conditional_prepared_eval(",
    source,
    fixed = TRUE
  )
  expect_gt(start, 0L)
  expect_gt(stop, start)
  evaluator <- substr(source, start, stop - 1L)

  expect_true(grepl(
    "if (np_conditional_density_prepared_context.degree_key_len > 0 &&\n      glp_degree == NULL)",
    evaluator,
    fixed = TRUE
  ))
  expect_true(grepl(
    "degree_work = np_conditional_density_prepared_context.eval_degree;",
    evaluator,
    fixed = TRUE
  ))
  expect_false(grepl(
    "rbw == NULL || glp_degree == NULL || out == NULL",
    evaluator,
    fixed = TRUE
  ))
})

test_that("scalar categorical-response CVML matches direct fixed and NN oracles", {
  fixture <- npcdens_y_side_fixture()
  eval_only <- getFromNamespace(".npcdensbw_eval_only", "np")
  cases <- list(
    fixed_unordered = list(
      y = fixture$y,
      bwtype = "fixed",
      bws = c(0.20, 0.22, 0.24, 0.10)
    ),
    generalized_nn_unordered = list(
      y = fixture$y,
      bwtype = "generalized_nn",
      bws = c(0.20, 12, 13, 0.10)
    )
  )

  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    bw <- npcdensbw(
      xdat = fixture$x,
      ydat = case$y,
      regtype = "lc",
      bwmethod = "cv.ml",
      bwtype = case$bwtype,
      bws = case$bws,
      bandwidth.compute = FALSE
    )
    native <- eval_only(fixture$x, case$y, bw)$objective
    oracle <- npcdens_y_side_oracle(fixture$x, case$y, bw)

    expect_true(is.finite(native), info = case_name)
    expect_equal(native, oracle, tolerance = 2e-12, info = case_name)
  }
})

test_that("prepared scalar categorical CVML retains the canonical finite penalty", {
  n <- 64L
  x <- data.frame(x = seq(0, 1, length.out = n))
  y <- data.frame(y = factor(rep(letters[1:5], length.out = n)))
  bw <- npcdensbw(
    xdat = x, ydat = y, regtype = "lc", bwmethod = "cv.ml",
    bwtype = "fixed", bws = c(0.4, 0.00025),
    bandwidth.compute = FALSE
  )
  oracle <- getFromNamespace(".npcdensbw_eval_only", "np")(
    xdat = x, ydat = y, bws = bw,
    invalid.penalty = "baseline", penalty.multiplier = 10
  )
  prep <- getFromNamespace(".npcdensbw_prepared_prepare_args", "np")(
    xdat = x, ydat = y, bws = bw, start.bw = c(bw$ybw, bw$xbw),
    invalid.penalty = "baseline", penalty.multiplier = 10,
    degree.search = TRUE
  )
  prepare <- getFromNamespace("npPreparedObjectivePrepareConditionalDensity", "np")
  evaluate <- getFromNamespace("npPreparedObjectiveEvalConditionalDensity", "np")
  destroy <- getFromNamespace("npPreparedObjectiveDestroyConditionalDensity", "np")

  prepared <- prepare(
    c.uno = prep$c.uno, c.ord = prep$c.ord, c.con = prep$c.con,
    u.uno = prep$u.uno, u.ord = prep$u.ord, u.con = prep$u.con,
    mysd = prep$mysd, myopti = prep$myopti, myoptd = prep$myoptd,
    rbw = prep$rbw, penalty.mode = prep$penalty_mode,
    penalty.multiplier = prep$penalty_multiplier, degree = prep$degree,
    bernstein = prep$bernstein, basis = prep$basis, regtype = prep$regtype,
    cxkerlb = prep$cxkerlb, cxkerub = prep$cxkerub,
    cykerlb = prep$cykerlb, cykerub = prep$cykerub
  )
  on.exit(if (isTRUE(prepared)) destroy(), add = TRUE)
  first <- evaluate(bw = as.double(prep$rbw), degree = as.integer(prep$degree))
  second <- evaluate(bw = as.double(prep$rbw), degree = as.integer(prep$degree))

  expect_true(prepared)
  expect_lt(oracle$objective, -1e3)
  expect_gt(oracle$objective, -1e7)
  expect_equal(first[[1L]], oracle$objective, tolerance = 2e-12)
  expect_identical(second, first)
})

test_that("Powell and native MADS initialize scalar categorical-response CVML", {
  fixture <- npcdens_y_side_fixture()
  x <- fixture$x[seq_len(32L), c("x1", "x2"), drop = FALSE]
  y <- fixture$y[seq_len(32L), , drop = FALSE]

  for (solver in c("powell", "mads")) {
    bw <- npcdensbw(
      xdat = x,
      ydat = y,
      regtype = "lc",
      bwmethod = "cv.ml",
      bwtype = "fixed",
      bwsolver = solver,
      nmulti = 1L
    )
    expect_true(is.finite(bw$fval), info = solver)
  }
})
