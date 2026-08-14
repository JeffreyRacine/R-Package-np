.npRmpi_bcast_cmd_expr <- getFromNamespace(".npRmpi_bcast_cmd_expr", "npRmpi")
.npRmpi_autodispatch_call <- getFromNamespace(".npRmpi_autodispatch_call", "npRmpi")
.npRmpi_manual_distributed_call <- getFromNamespace(".npRmpi_manual_distributed_call", "npRmpi")
.npRmpi_bcast_robj_by_name <- getFromNamespace(".npRmpi_bcast_robj_by_name", "npRmpi")
.npRmpi_eval_without_dispatch <- getFromNamespace(".npRmpi_eval_without_dispatch", "npRmpi")
.npRmpi_autodispatch_eval_arg <- getFromNamespace(".npRmpi_autodispatch_eval_arg", "npRmpi")
.npRmpi_autodispatch_cleanup <- getFromNamespace(".npRmpi_autodispatch_cleanup", "npRmpi")
.npRmpi_distributed_call_impl <- getFromNamespace(".npRmpi_distributed_call_impl", "npRmpi")
.npRmpi_autodispatch_materialize_call <- getFromNamespace(".npRmpi_autodispatch_materialize_call", "npRmpi")
.npRmpi_bootstrap_compute_payload <- getFromNamespace(".npRmpi_bootstrap_compute_payload", "npRmpi")
.npRmpi_rm_existing <- getFromNamespace(".npRmpi_rm_existing", "npRmpi")
.np_eval_bws_call_arg <- getFromNamespace(".np_eval_bws_call_arg", "npRmpi")
.npRmpi_autodispatch_target_args <- getFromNamespace(".npRmpi_autodispatch_target_args", "npRmpi")
.npRmpi_autodispatch_replace_tmps <- getFromNamespace(".npRmpi_autodispatch_replace_tmps", "npRmpi")
.npRmpi_autodispatch_sanitize_object <- getFromNamespace(".npRmpi_autodispatch_sanitize_object", "npRmpi")
.npRmpi_autodispatch_remote_ref <- getFromNamespace(".npRmpi_autodispatch_remote_ref", "npRmpi")
.npRmpi_autodispatch_fingerprint <- getFromNamespace(".npRmpi_autodispatch_fingerprint", "npRmpi")
.npRmpi_autodispatch_raw_md5 <- getFromNamespace(".npRmpi_autodispatch_raw_md5", "npRmpi")
.npRmpi_autodispatch_tag_result <- getFromNamespace(".npRmpi_autodispatch_tag_result", "npRmpi")
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

  expect_match(cleanup.body, "\\.npRmpi_lease_run_lifecycle")
  expect_match(cleanup.body, "autodispatch.lifecycle.retire", fixed = TRUE)
  expect_match(cleanup.body, "\\.npRmpi_rm_existing\\(ordinary, envir = \\.GlobalEnv\\)")
  expect_match(cleanup.body, "get\\(\"\\.npRmpi_rm_existing\"")
  expect_match(cleanup.body, "asNamespace\\(\"npRmpi\"\\)")
  expect_match(cleanup.body, "TMPS, envir = \\.GlobalEnv")
  expect_match(impl.body, "\\.npRmpi_autodispatch_cleanup\\(prepared\\$tmpnames, comm = comm\\)")
  expect_match(boot.body, "\\.npRmpi_rm_existing\\(tmp, envir = \\.GlobalEnv\\)")
  expect_match(boot.body, "get\\(\"\\.npRmpi_rm_existing\"")
  expect_match(boot.body, "asNamespace\\(\"npRmpi\"\\)")
  expect_match(boot.body, "TMP, envir = \\.GlobalEnv")
})

test_that("autodispatch return rewriting covers all coordinator replacements", {
  impl.body <- paste(deparse(body(.npRmpi_distributed_call_impl), width.cutoff = 500L), collapse = " ")
  sanitize.body <- paste(deparse(body(.npRmpi_autodispatch_sanitize_object), width.cutoff = 500L), collapse = " ")
  eval.body <- paste(deparse(body(getFromNamespace(".npRmpi_spmd_eval_payload", "npRmpi")), width.cutoff = 500L), collapse = " ")

  expect_match(impl.body, "prepared\\$lease.replacements")
  expect_match(impl.body, "prepublish\\.names = names\\(prepared\\$prepublish\\)")
  expect_match(eval.body, "prepublish\\.names <- payload\\$prepublish\\.names")
  expect_match(impl.body, "\\.npRmpi_autodispatch_sanitize_object\\(result, tmpvals = tmpreplace\\)")
  expect_match(impl.body, "\\.npRmpi_autodispatch_replace_tmps\\(result, tmpvals = tmpreplace\\)")
  expect_match(sanitize.body, "\\.npRmpi_autodispatch_replace_tmp_calls\\(x, tmpvals = tmpvals\\)")
})

test_that("autodispatch prepublishes large implicit formula data", {
  withr::local_options(npRmpi.autodispatch.arg.broadcast.threshold.regression = 1L)

  env <- new.env(parent = .GlobalEnv)
  env$x <- seq_len(10)
  env$y <- env$x + 1

  prepared <- .npRmpi_autodispatch_materialize_call(
    quote(npregbw(y ~ x, regtype = "ll")),
    caller_env = env
  )

  expect_true(any(grepl("^\\.__npRmpi_autod_data_", prepared$tmpnames)))
  expect_true(any(grepl("^\\.__npRmpi_autod_data_", names(prepared$prepublish))))
  expect_false(any(grepl("^\\.__npRmpi_autod_data_", names(prepared$tmpvals))))
})

test_that("autodispatch reuses semiparametric remote bandwidth references", {
  withr::local_options(npRmpi.autodispatch.arg.broadcast.threshold.regression = 1L)

  reset <- getFromNamespace(".npRmpi_lease_reset_local", "npRmpi")
  plan <- getFromNamespace(".npRmpi_lease_publication_plan", "npRmpi")
  prepare <- getFromNamespace(".npRmpi_lease_prepare_local", "npRmpi")
  commit <- getFromNamespace(".npRmpi_lease_commit_local", "npRmpi")
  reset()
  on.exit(reset(), add = TRUE)

  env <- new.env(parent = .GlobalEnv)
  env$sibw <- structure(
    list(call = quote(npindexbw(y ~ x1 + x2, data = mydat)),
         formula = y ~ x1 + x2,
         ballast = seq_len(100)),
    class = "sibandwidth"
  )
  si.plan <- plan(quote(npindexbw(xdat = x, ydat = y)))
  prepare(env$sibw, si.plan)
  commit(si.plan)
  env$sibw <- .npRmpi_autodispatch_tag_result(env$sibw, publication = si.plan)
  env$scbw <- structure(
    list(call = quote(npscoefbw(y ~ x | z, data = mydat)),
         formula = y ~ x | z,
         ballast = seq_len(100)),
    class = "scbandwidth"
  )
  sc.plan <- plan(quote(npscoefbw(xdat = x, ydat = y, zdat = z)))
  prepare(env$scbw, sc.plan)
  commit(sc.plan)
  env$scbw <- .npRmpi_autodispatch_tag_result(env$scbw, publication = sc.plan)

  si <- .npRmpi_autodispatch_materialize_call(
    quote(npindex(bws = sibw, gradients = FALSE)),
    caller_env = env
  )
  sc <- .npRmpi_autodispatch_materialize_call(
    quote(npscoef(bws = scbw, gradients = FALSE)),
    caller_env = env
  )

  expect_length(si$lease.bindings, 1L)
  expect_length(sc$lease.bindings, 1L)
  expect_identical(unname(si$lease.bindings[[1L]]), si.plan$id)
  expect_identical(unname(sc$lease.bindings[[1L]]), sc.plan$id)
  expect_identical(as.character(si$call$bws), names(si$lease.bindings))
  expect_identical(as.character(sc$call$bws), names(sc$lease.bindings))
  expect_identical(names(si$lease.replacements), names(si$lease.bindings))
  expect_identical(names(sc$lease.replacements), names(sc$lease.bindings))
  expect_identical(si$lease.replacements[[1L]],
                   getFromNamespace(".npRmpi_autodispatch_untag", "npRmpi")(env$sibw))
  expect_identical(sc$lease.replacements[[1L]],
                   getFromNamespace(".npRmpi_autodispatch_untag", "npRmpi")(env$scbw))
  expect_false(names(si$lease.replacements) %in% names(si$tmpvals))
  expect_false(names(sc$lease.replacements) %in% names(sc$tmpvals))
  expect_false(names(si$lease.replacements) %in% names(si$prepublish))
  expect_false(names(sc$lease.replacements) %in% names(sc$prepublish))
  expect_false(any(vapply(si$prepublish, inherits, logical(1), "sibandwidth")))
  expect_false(any(vapply(sc$prepublish, inherits, logical(1), "scbandwidth")))
})

test_that("autodispatch remote bandwidth fingerprints use an MD5 digest", {
  bw <- structure(
    list(call = quote(npindexbw(y ~ x1 + x2, data = mydat)),
         formula = y ~ x1 + x2,
         bw = c(0.4, 0.6),
         ballast = seq_len(100)),
    class = "sibandwidth"
  )

  fp <- .npRmpi_autodispatch_fingerprint(bw)
  expect_match(fp, "^md5:[0-9]+:[0-9a-f]{32}$")

  bw_changed <- bw
  bw_changed$bw[1L] <- bw_changed$bw[1L] + 0.01
  expect_false(identical(fp, .npRmpi_autodispatch_fingerprint(bw_changed)))

  body_text <- paste(
    deparse(body(.npRmpi_autodispatch_raw_md5), width.cutoff = 500L),
    collapse = " "
  )
  expect_true(grepl("tools::md5sum", body_text, fixed = TRUE))
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

test_that("autodispatch sanitation preserves explicit NULL arguments and later fields", {
  call_in <- quote(npudistbw(
    dat = .__npRmpi_autod_dat_1,
    gdat = NULL,
    bws = .__npRmpi_autod_bws_3,
    bandwidth.compute = FALSE,
    do.full.integral = TRUE
  ))
  tmpvals <- list(
    .__npRmpi_autod_dat_1 = data.frame(x = 1:3),
    .__npRmpi_autod_bws_3 = 0.25
  )

  sanitized <- .npRmpi_autodispatch_sanitize_object(
    list(call = call_in, marker = "retained"),
    tmpvals = tmpvals
  )
  call.args <- as.list(sanitized$call)

  expect_identical(names(call.args), names(as.list(call_in)))
  expect_identical(call.args$dat, tmpvals$.__npRmpi_autod_dat_1)
  expect_true("gdat" %in% names(call.args))
  expect_null(call.args$gdat)
  expect_identical(call.args$bws, 0.25)
  expect_identical(call.args$bandwidth.compute, FALSE)
  expect_identical(call.args$do.full.integral, TRUE)
  expect_identical(sanitized$marker, "retained")

  pair_in <- as.pairlist(as.list(call_in)[-1L])
  pair_out <- .npRmpi_autodispatch_replace_tmps(pair_in, tmpvals = tmpvals)
  pair.args <- as.list(pair_out)
  expect_identical(names(pair.args), names(as.list(pair_in)))
  expect_true("gdat" %in% names(pair.args))
  expect_null(pair.args$gdat)
  expect_identical(pair.args$do.full.integral, TRUE)
})

test_that("npudistbw autodispatch preserves explicit NULL training-grid calls", {
  if (!spawn_mpi_slaves(1L))
    skip("Could not spawn MPI slave")
  on.exit(close_mpi_slaves(), add = TRUE)

  old.options <- options(np.messages = FALSE)
  on.exit(options(old.options), add = TRUE)

  set.seed(2026073011L)
  values <- seq_len(8L)
  dat <- data.frame(
    o1 = ordered(sample(values, 480L, replace = TRUE), levels = values)
  )
  bw <- npudistbw(
    dat = dat,
    gdat = NULL,
    bws = 0.19,
    bwmethod = "cv.cdf",
    bwtype = "fixed",
    okertype = "wangvanryzin",
    do.full.integral = TRUE,
    bandwidth.compute = FALSE
  )
  objective <- npudistbw(
    dat = dat,
    gdat = NULL,
    bws = bw,
    do.full.integral = TRUE,
    eval.only = TRUE
  )

  expect_s3_class(bw, "dbandwidth")
  expect_s3_class(objective, "dbandwidth")
  expect_true(is.finite(objective$fval))
  expect_identical(objective$bw, bw$bw)
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
