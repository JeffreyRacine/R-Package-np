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
  local_eval <- getFromNamespace(".npRmpi_with_local_regression", "npRmpi")
  xweights <- local_eval(npksum(
    txdat = x, bws = bw$xbw, bwtype = bw$type,
    ckertype = bw$cxkertype, ckerorder = bw$cxkerorder,
    ukertype = bw$uxkertype, okertype = bw$oxkertype,
    return.kernel.weights = TRUE
  ))$kw
  yweights <- local_eval(npksum(
    txdat = y, bws = bw$ybw, bwtype = bw$type,
    ukertype = bw$uykertype, okertype = bw$oykertype,
    return.kernel.weights = TRUE
  ))$kw
  diag(xweights) <- 0
  diag(yweights) <- 0
  sum(log(colSums(xweights * yweights) / colSums(xweights)))
}

test_that("conditional CVML Y-side ownership covers every searchable stream", {
  source <- npcdens_y_side_source("np.c")
  skip_if(is.null(source), "source file src/np.c unavailable")

  expect_true(grepl(
    "static int np_conditional_density_objective_needs_y_side(",
    source, fixed = TRUE
  ))
  expect_match(
    source,
    "\\(degree_search \\|\\|\\s*np_conditional_density_cvml_stream_engine_supported\\(\\)\\)\\);",
    perl = TRUE
  )
  expect_equal(
    lengths(regmatches(
      source,
      gregexpr(
        "need_y_side = np_conditional_density_objective_needs_y_side(\n    ibwmfunc, degree_search);",
        source, fixed = TRUE
      )
    )),
    1L
  )
  expect_false(grepl(
    "need_y_side = (ibwmfunc == CBWM_CVLS) || ((ibwmfunc == CBWM_CVML) && (np_lp_engine_extern == NP_LP_ENGINE_GENERAL));",
    source, fixed = TRUE
  ))
})

test_that("categorical kernel rows reject missing category metadata", {
  source <- npcdens_y_side_source("jksum.c")
  skip_if(is.null(source), "source file src/jksum.c unavailable")

  start <- regexpr(
    "static int NP_OUTER_PACK_ADJACENT_HOT_ALIGN kernel_weighted_sum_np_ctx_ex(",
    source, fixed = TRUE
  )
  stop <- regexpr("int kernel_weighted_sum_np_ctx(", source,
                  fixed = TRUE)
  expect_gt(start, 0L)
  expect_gt(stop, start)
  owner <- substr(source, start, stop - 1L)
  expect_true(grepl(
    "((num_reg_unordered > 0) || (num_reg_ordered > 0))",
    owner, fixed = TRUE
  ))
  expect_true(grepl("(num_categories == NULL)", owner, fixed = TRUE))
})

test_that("native CVML degree buffers use bounded prepared scratch", {
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
  expect_true(grepl(
    "if (np_conditional_density_prepared_context.degree_key_len > 0)\n      MPI_Bcast(degree_work,",
    evaluator,
    fixed = TRUE
  ))
  expect_true(grepl(
    paste0(
      "if (np_conditional_density_prepared_context.degree_search ||\n",
      "      degree_refresh_needed)\n",
      "    degree_refresh_ok = np_mpi_comm1_all_ok(degree_refresh_ok);"
    ),
    evaluator,
    fixed = TRUE
  ))
  expect_false(grepl(
    "rbw == NULL || glp_degree == NULL || out == NULL",
    evaluator,
    fixed = TRUE
  ))
})

test_that("scalar categorical-response CVML matches the direct fixed oracle", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old_options <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old_options), add = TRUE)

  fixture <- npcdens_y_side_fixture()
  eval_only <- getFromNamespace(".npcdensbw_eval_only", "npRmpi")
  cases <- list(
    fixed = list(bwtype = "fixed", bws = c(0.20, 0.22, 0.24, 0.10))
  )

  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    bw <- npcdensbw(
      xdat = fixture$x,
      ydat = fixture$y,
      regtype = "lc",
      bwmethod = "cv.ml",
      bwtype = case$bwtype,
      bws = case$bws,
      bandwidth.compute = FALSE
    )
    native <- eval_only(fixture$x, fixture$y, bw)$objective
    oracle <- npcdens_y_side_oracle(fixture$x, fixture$y, bw)

    expect_true(is.finite(native), info = case_name)
    expect_equal(native, oracle, tolerance = 2e-12, info = case_name)
  }
})

test_that("prepared scalar categorical CVML retains the canonical finite penalty", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old_options <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old_options), add = TRUE)

  n <- 64L
  x <- data.frame(x = seq(0, 1, length.out = n))
  y <- data.frame(y = factor(rep(letters[1:5], length.out = n)))
  local_eval <- getFromNamespace(".npRmpi_with_local_regression", "npRmpi")
  bw <- local_eval(npcdensbw(
    xdat = x, ydat = y, regtype = "lc", bwmethod = "cv.ml",
    bwtype = "fixed", bws = c(0.4, 0.00025),
    bandwidth.compute = FALSE
  ))
  oracle <- getFromNamespace(".npcdensbw_eval_only", "npRmpi")(
    xdat = x, ydat = y, bws = bw,
    invalid.penalty = "baseline", penalty.multiplier = 10
  )
  prep <- getFromNamespace(".npcdensbw_prepared_prepare_args", "npRmpi")(
    xdat = x, ydat = y, bws = bw, start.bw = c(bw$ybw, bw$xbw),
    invalid.penalty = "baseline", penalty.multiplier = 10,
    degree.search = TRUE
  )

  cmd <- substitute(local({
    ns <- asNamespace("npRmpi")
    prepare <- get(
      "npRmpiPreparedObjectivePrepareConditionalDensity", envir = ns
    )
    evaluate <- get(
      "npRmpiPreparedObjectiveEvalConditionalDensity", envir = ns
    )
    destroy <- get(
      "npRmpiPreparedObjectiveDestroyConditionalDensity", envir = ns
    )
    p <- PREP
    prepared <- prepare(
      c.uno = p$c.uno, c.ord = p$c.ord, c.con = p$c.con,
      u.uno = p$u.uno, u.ord = p$u.ord, u.con = p$u.con,
      mysd = p$mysd, myopti = p$myopti, myoptd = p$myoptd,
      rbw = p$rbw, penalty.mode = p$penalty_mode,
      penalty.multiplier = p$penalty_multiplier, degree = p$degree,
      bernstein = p$bernstein, basis = p$basis, regtype = p$regtype,
      cxkerlb = p$cxkerlb, cxkerub = p$cxkerub,
      cykerlb = p$cykerlb, cykerub = p$cykerub
    )
    mpi.barrier(1L)
    first <- evaluate(bw = as.double(p$rbw), degree = as.integer(p$degree))
    second <- evaluate(bw = as.double(p$rbw), degree = as.integer(p$degree))
    mpi.barrier(1L)
    destroy()
    list(prepared = prepared, first = first, second = second)
  }), list(PREP = prep))
  distributed <- getFromNamespace(".npRmpi_bcast_cmd_expr", "npRmpi")(
    cmd, comm = 1L, caller.execute = TRUE
  )

  expect_true(distributed$prepared)
  expect_lt(oracle$objective, -1e3)
  expect_gt(oracle$objective, -1e7)
  expect_equal(distributed$first[[1L]], oracle$objective, tolerance = 2e-12)
  expect_identical(distributed$second, distributed$first)
})
