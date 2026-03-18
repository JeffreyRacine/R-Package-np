.npRmpi_bcast_cmd_expr <- getFromNamespace(".npRmpi_bcast_cmd_expr", "npRmpi")
.npRmpi_autodispatch_call <- getFromNamespace(".npRmpi_autodispatch_call", "npRmpi")
.npRmpi_manual_distributed_call <- getFromNamespace(".npRmpi_manual_distributed_call", "npRmpi")
.npRmpi_bcast_robj_by_name <- getFromNamespace(".npRmpi_bcast_robj_by_name", "npRmpi")
.npRmpi_eval_without_dispatch <- getFromNamespace(".npRmpi_eval_without_dispatch", "npRmpi")
.npRmpi_autodispatch_eval_arg <- getFromNamespace(".npRmpi_autodispatch_eval_arg", "npRmpi")
.npRmpi_autodispatch_cleanup <- getFromNamespace(".npRmpi_autodispatch_cleanup", "npRmpi")
.npRmpi_distributed_call_impl <- getFromNamespace(".npRmpi_distributed_call_impl", "npRmpi")
.npRmpi_bootstrap_compute_payload <- getFromNamespace(".npRmpi_bootstrap_compute_payload", "npRmpi")
.npRmpi_rm_existing <- getFromNamespace(".npRmpi_rm_existing", "npRmpi")
.np_eval_bws_call_arg <- getFromNamespace(".np_eval_bws_call_arg", "npRmpi")
.npRmpi_autodispatch_target_args <- getFromNamespace(".npRmpi_autodispatch_target_args", "npRmpi")
.npRmpi_autodispatch_replace_tmps <- getFromNamespace(".npRmpi_autodispatch_replace_tmps", "npRmpi")
.npRmpi_is_missing_call_arg <- getFromNamespace(".npRmpi_is_missing_call_arg", "npRmpi")

test_that(".npRmpi_bcast_cmd_expr forwards command expression structurally", {
  env <- new.env(parent = environment())
  env$seen <- NULL
  env$mpi.bcast.cmd <- function(cmd, comm = 1L, caller.execute = TRUE) {
    env$seen <- list(expr = substitute(cmd),
                     value = cmd,
                     comm = comm,
                     caller.execute = caller.execute)
    "OK"
  }

  out <- evalq(.npRmpi_bcast_cmd_expr(quote(x <- 1L), comm = 3L, caller.execute = FALSE), envir = env)

  expect_identical(out, "OK")
  expect_true(is.list(env$seen))
  expect_true(is.call(env$seen$expr))
  expect_identical(env$seen$expr, quote(x <- 1L))
  expect_identical(env$seen$value, 1L)
  expect_identical(env$seen$comm, 3L)
  expect_identical(env$seen$caller.execute, FALSE)
})

test_that(".npRmpi_bcast_cmd_expr resolves mpi.bcast.cmd from caller frame", {
  fn.body <- paste(deparse(body(.npRmpi_bcast_cmd_expr), width.cutoff = 500L), collapse = " ")
  expect_match(fn.body, "get\\(\"mpi\\.bcast\\.cmd\", envir = parent\\.frame\\(\\), mode = \"function\", inherits = TRUE\\)")
})

test_that(".npRmpi_autodispatch_call delegates to shared distributed-call helper", {
  fn.body <- paste(deparse(body(.npRmpi_autodispatch_call), width.cutoff = 500L), collapse = " ")
  expect_match(fn.body, "\\.npRmpi_distributed_call_impl\\(mc = mc, caller_env = caller_env, comm = comm, warn_nested = TRUE\\)")
})

test_that(".npRmpi_manual_distributed_call delegates to shared distributed-call helper", {
  fn.body <- paste(deparse(body(.npRmpi_manual_distributed_call), width.cutoff = 500L), collapse = " ")
  expect_match(fn.body, "\\.npRmpi_distributed_call_impl\\(mc = mc, caller_env = caller_env, comm = comm, warn_nested = FALSE\\)")
})

test_that("autodispatch eval helpers route through shared command executor", {
  bcast.body <- paste(deparse(body(.npRmpi_bcast_robj_by_name), width.cutoff = 500L), collapse = " ")
  nodisp.body <- paste(deparse(body(.npRmpi_eval_without_dispatch), width.cutoff = 500L), collapse = " ")
  arg.body <- paste(deparse(body(.npRmpi_autodispatch_eval_arg), width.cutoff = 500L), collapse = " ")

  expect_match(bcast.body, "\\.npRmpi_eval_scmd\\(call, envir = caller_env\\)")
  expect_match(nodisp.body, "\\.npRmpi_eval_scmd\\(mc\\.eval, envir = caller_env\\)")
  expect_match(arg.body, "\\.npRmpi_eval_scmd\\(expr, envir = caller_env\\)")
})

test_that("autodispatch uses safe cleanup helper for temporary symbols", {
  cleanup.body <- paste(deparse(body(.npRmpi_autodispatch_cleanup), width.cutoff = 500L), collapse = " ")
  impl.body <- paste(deparse(body(.npRmpi_distributed_call_impl), width.cutoff = 500L), collapse = " ")
  boot.body <- paste(deparse(body(.npRmpi_bootstrap_compute_payload), width.cutoff = 500L), collapse = " ")

  expect_match(cleanup.body, "\\.npRmpi_rm_existing\\(tmpnames, envir = \\.GlobalEnv\\)")
  expect_match(cleanup.body, "get\\(\"\\.npRmpi_rm_existing\"")
  expect_match(cleanup.body, "asNamespace\\(\"npRmpi\"\\)")
  expect_match(cleanup.body, "TMPS, envir = \\.GlobalEnv")
  expect_match(impl.body, "\\.npRmpi_autodispatch_cleanup\\(prepared\\$tmpnames, comm = comm\\)")
  expect_match(boot.body, "\\.npRmpi_rm_existing\\(tmp, envir = \\.GlobalEnv\\)")
  expect_match(boot.body, "get\\(\"\\.npRmpi_rm_existing\"")
  expect_match(boot.body, "asNamespace\\(\"npRmpi\"\\)")
  expect_match(boot.body, "TMP, envir = \\.GlobalEnv")
})

test_that("autodispatch return rewriting covers prepublished temporary arguments", {
  impl.body <- paste(deparse(body(.npRmpi_distributed_call_impl), width.cutoff = 500L), collapse = " ")
  sanitize.body <- paste(deparse(body(.npRmpi_autodispatch_sanitize_object), width.cutoff = 500L), collapse = " ")

  expect_match(impl.body, "tmpreplace <- c\\(prepared\\$tmpvals, prepared\\$prepublish\\)")
  expect_match(impl.body, "\\.npRmpi_autodispatch_sanitize_object\\(result, tmpvals = tmpreplace\\)")
  expect_match(impl.body, "\\.npRmpi_autodispatch_replace_tmps\\(result, tmpvals = tmpreplace\\)")
  expect_match(sanitize.body, "\\.npRmpi_autodispatch_replace_tmp_calls\\(x, tmpvals = tmpvals\\)")
})

test_that("autodispatch target argument set covers gdat alias", {
  args <- .npRmpi_autodispatch_target_args()
  expect_true(is.character(args))
  expect_true("gdata" %in% args)
  expect_true("gdat" %in% args)
})

test_that("autodispatch tmp replacement handles calls with missing arguments", {
  call_in <- quote(npplregbw(formula = y ~ x, data = , bws = .__npRmpi_autod_bws_1))
  call_out <- .npRmpi_autodispatch_replace_tmps(
    call_in,
    tmpvals = list(".__npRmpi_autod_bws_1" = 7L)
  )

  out_list <- as.list(call_out)
  expect_true(is.call(call_out))
  expect_identical(out_list[[1L]], as.name("npplregbw"))
  expect_true(.npRmpi_is_missing_call_arg(out_list[[3L]]))
  expect_identical(out_list[[4L]], 7L)
})

test_that("npudist(bws=...) resolves large autodispatch temporary call arguments", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  withr::local_options(npRmpi.autodispatch.arg.broadcast.threshold = 1L)

  data("faithful")
  bw <- npudistbw(dat = faithful, bws = c(0.5, 5), bandwidth.compute = FALSE)

  resolved_dat <- .np_eval_bws_call_arg(bw, "dat")
  expect_true(is.data.frame(resolved_dat))
  expect_equal(nrow(resolved_dat), nrow(faithful))

  fit <- npudist(bws = bw)
  expect_s3_class(fit, "npdistribution")
  expect_equal(length(fitted(fit)), nrow(faithful))
})

test_that(".npRmpi_rm_existing removes only existing names", {
  env <- new.env(parent = emptyenv())
  env$foo <- 1L
  expect_silent(.npRmpi_rm_existing(c("foo", "bar"), envir = env))
  expect_false(exists("foo", envir = env, inherits = FALSE))
})
