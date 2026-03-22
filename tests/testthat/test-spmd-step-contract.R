run_spmd_subprocess <- function(lines, timeout = 60L, env = character()) {
  script <- tempfile("npRmpi-spmd-", fileext = ".R")
  writeLines(lines, script, useBytes = TRUE)
  on.exit(unlink(script), add = TRUE)

  cmd <- file.path(R.home("bin"), "Rscript")
  out <- suppressWarnings(system2(
    cmd,
    c("--no-save", script),
    stdout = TRUE,
    stderr = TRUE,
    timeout = timeout,
    env = env
  ))
  status <- attr(out, "status")
  if (is.null(status))
    status <- 0L
  list(status = as.integer(status), output = out)
}

test_that("SPMD envelope sequence contract fails fast on mismatch", {
  make.env <- getFromNamespace(".npRmpi_spmd_make_envelope", "npRmpi")
  exec.env <- getFromNamespace(".npRmpi_spmd_execute_local", "npRmpi")
  try.exec <- getFromNamespace(".npRmpi_spmd_try_execute_local", "npRmpi")
  set.seq <- getFromNamespace(".npRmpi_spmd_seq_set", "npRmpi")

  old.seq <- getOption("npRmpi.spmd.seq_id", 0L)
  on.exit(options(npRmpi.spmd.seq_id = old.seq), add = TRUE)

  set.seq(0L)
  env1 <- make.env(opcode = "spmd.ping", timeout_class = "smoke")
  ok1 <- exec.env(env1, payload = list(label = "first"))
  expect_true(is.list(ok1) && isTRUE(ok1$ok))
  expect_identical(ok1$ack$status, "ACK")
  expect_identical(ok1$ack$seq_id, 1L)

  env.bad <- make.env(opcode = "spmd.ping", seq_id = 3L, timeout_class = "smoke")
  expect_error(
    exec.env(env.bad, payload = list(label = "bad")),
    "sequence mismatch"
  )

  err <- try.exec(env.bad, payload = list(label = "bad"))
  expect_true(is.list(err) && !isTRUE(err$ok))
  expect_identical(err$ack$status, "ERR")
  expect_match(err$error, "sequence mismatch")
})

test_that("SPMD tiny smoke opcode runs in session mode subprocess", {
  skip_on_cran()
  res <- run_spmd_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE)",
      "f <- getFromNamespace('.npRmpi_spmd_tiny_smoke', 'npRmpi')",
      "ans <- f(label='session-smoke')",
      "stopifnot(is.list(ans), isTRUE(ans$ok), identical(ans$ack$status, 'ACK'))",
      "cat('SPMD_TINY_SESSION_OK\\n')"
    ),
    timeout = 90L
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SPMD_TINY_SESSION_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("SPMD step execution uses collective ACK in session subprocess", {
  skip_on_cran()
  res <- run_spmd_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE)",
      "mk <- getFromNamespace('.npRmpi_spmd_make_envelope', 'npRmpi')",
      "run <- getFromNamespace('.npRmpi_spmd_execute_step', 'npRmpi')",
      "env <- mk(opcode='spmd.ping', timeout_class='smoke')",
      "ans <- run(envelope=env, payload=list(label='collective-ack'), comm=1L, where='spmd-collective')",
      "stopifnot(is.list(ans), isTRUE(ans$ok), identical(ans$ack$status, 'ACK'))",
      "cat('SPMD_COLLECTIVE_ACK_OK\\n')"
    ),
    timeout = 90L
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SPMD_COLLECTIVE_ACK_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("SPMD step diverged sequence fails fast with ACK mismatch diagnostics", {
  skip_on_cran()
  res <- run_spmd_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE)",
      "mk <- getFromNamespace('.npRmpi_spmd_make_envelope', 'npRmpi')",
      "run <- getFromNamespace('.npRmpi_spmd_execute_step', 'npRmpi')",
      "mpi.bcast.cmd(if (mpi.comm.rank(1L) == 1L) options(npRmpi.spmd.seq_id = 5L), caller.execute=TRUE)",
      "env <- mk(opcode='spmd.ping', timeout_class='smoke')",
      "err <- try(run(envelope=env, payload=list(label='diverge'), comm=1L, where='spmd-diverge'), silent=TRUE)",
      "if (!inherits(err, 'try-error')) stop('expected divergence failure')",
      "msg <- as.character(err)",
      "if (!any(grepl('ACK mismatch', msg, fixed=TRUE))) stop('missing ACK mismatch diagnostic')",
      "cat('SPMD_DIVERGENCE_FAILFAST_OK\\n')"
    ),
    timeout = 90L
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SPMD_DIVERGENCE_FAILFAST_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("SPMD opcode selection tags LL/LP CV routes for core bw families", {
  opcode.fun <- getFromNamespace(".npRmpi_spmd_opcode_from_call", "npRmpi")

  mc.reg <- quote(npregbw(xdat = x, ydat = y, regtype = "ll", bwmethod = "cv.ls"))
  op.reg <- opcode.fun(mc = mc.reg, caller_env = environment())
  expect_identical(op.reg, "autodispatch.npregbw.cv_lllp")

  mc.sc <- quote(npscoefbw(xdat = x, ydat = y, zdat = z, regtype = "lp", degree = 1L, bwmethod = "cv.aic"))
  op.sc <- opcode.fun(mc = mc.sc, caller_env = environment())
  expect_identical(op.sc, "autodispatch.npscoefbw.cv_lllp")

  mc.pl <- quote(npplregbw(xdat = x, ydat = y, zdat = z, regtype = "ll", bwmethod = "cv.ls"))
  op.pl <- opcode.fun(mc = mc.pl, caller_env = environment())
  expect_identical(op.pl, "autodispatch.npplregbw.cv_lllp")

  mc.si <- quote(npindexbw(xdat = x, ydat = y))
  op.si <- opcode.fun(mc = mc.si, caller_env = environment())
  expect_identical(op.si, "autodispatch.npindexbw.core")

  mc.si.est <- quote(npindex(bws = bw, gradients = FALSE))
  op.si.est <- opcode.fun(mc = mc.si.est, caller_env = environment())
  expect_identical(op.si.est, "autodispatch.npindex.core")

  mc.reg.est <- quote(npreg(bws = bw, gradients = FALSE))
  op.reg.est <- opcode.fun(mc = mc.reg.est, caller_env = environment())
  expect_identical(op.reg.est, "autodispatch.npreg.core")

  mc.sc.est <- quote(npscoef(bws = bw, gradients = FALSE))
  op.sc.est <- opcode.fun(mc = mc.sc.est, caller_env = environment())
  expect_identical(op.sc.est, "autodispatch.npscoef.core")

  mc.pl.est <- quote(npplreg(bws = bw, gradients = FALSE))
  op.pl.est <- opcode.fun(mc = mc.pl.est, caller_env = environment())
  expect_identical(op.pl.est, "autodispatch.npplreg.core")

  mc.ud.est <- quote(npudens(tdat = x, bws = bw))
  op.ud.est <- opcode.fun(mc = mc.ud.est, caller_env = environment())
  expect_identical(op.ud.est, "autodispatch.npudens.core")

  mc.ui.est <- quote(npudist(tdat = x, bws = bw))
  op.ui.est <- opcode.fun(mc = mc.ui.est, caller_env = environment())
  expect_identical(op.ui.est, "autodispatch.npudist.core")

  mc.cd.est <- quote(npcdens(txdat = x, tydat = y, bws = bw))
  op.cd.est <- opcode.fun(mc = mc.cd.est, caller_env = environment())
  expect_identical(op.cd.est, "autodispatch.npcdens.core")

  mc.ci.est <- quote(npcdist(txdat = x, tydat = y, bws = bw))
  op.ci.est <- opcode.fun(mc = mc.ci.est, caller_env = environment())
  expect_identical(op.ci.est, "autodispatch.npcdist.core")

  mc.qreg <- quote(npqreg(bws = bw, tau = 0.5))
  op.qreg <- opcode.fun(mc = mc.qreg, caller_env = environment())
  expect_identical(op.qreg, "autodispatch.npqreg.core")

  mc.conmode <- quote(npconmode(bws = bw))
  op.conmode <- opcode.fun(mc = mc.conmode, caller_env = environment())
  expect_identical(op.conmode, "autodispatch.npconmode.core")

  mc.ksum <- quote(npksum(txdat = x, bws = 0.5))
  op.ksum <- opcode.fun(mc = mc.ksum, caller_env = environment())
  expect_identical(op.ksum, "autodispatch.npksum.core")

  mc.iv <- quote(npregiv(y = y, z = z, w = w))
  op.iv <- opcode.fun(mc = mc.iv, caller_env = environment())
  expect_identical(op.iv, "autodispatch.npregiv.core")

  mc.ivd <- quote(npregivderiv(y = y, z = z, w = w))
  op.ivd <- opcode.fun(mc = mc.ivd, caller_env = environment())
  expect_identical(op.ivd, "autodispatch.npregivderiv.core")

  mc.cmt <- quote(npcmstest(model = m, xdat = x, ydat = y))
  op.cmt <- opcode.fun(mc = mc.cmt, caller_env = environment())
  expect_identical(op.cmt, "autodispatch.npcmstest.core")

  mc.qcmt <- quote(npqcmstest(model = m, xdat = x, ydat = y))
  op.qcmt <- opcode.fun(mc = mc.qcmt, caller_env = environment())
  expect_identical(op.qcmt, "autodispatch.npqcmstest.core")

  mc.deq <- quote(npdeneqtest(x = xdf, y = ydf))
  op.deq <- opcode.fun(mc = mc.deq, caller_env = environment())
  expect_identical(op.deq, "autodispatch.npdeneqtest.core")

  mc.dep <- quote(npdeptest(data.x = x, data.y = y))
  op.dep <- opcode.fun(mc = mc.dep, caller_env = environment())
  expect_identical(op.dep, "autodispatch.npdeptest.core")

  mc.sdep <- quote(npsdeptest(data = y))
  op.sdep <- opcode.fun(mc = mc.sdep, caller_env = environment())
  expect_identical(op.sdep, "autodispatch.npsdeptest.core")

  mc.sym <- quote(npsymtest(data = y))
  op.sym <- opcode.fun(mc = mc.sym, caller_env = environment())
  expect_identical(op.sym, "autodispatch.npsymtest.core")

  mc.uni <- quote(npunitest(data.x = x, data.y = y))
  op.uni <- opcode.fun(mc = mc.uni, caller_env = environment())
  expect_identical(op.uni, "autodispatch.npunitest.core")

  mc.noncv <- quote(npregbw(xdat = x, ydat = y, regtype = "lc", bwmethod = "cv.ls"))
  op.noncv <- opcode.fun(mc = mc.noncv, caller_env = environment())
  expect_identical(op.noncv, "autodispatch.npregbw")
})

test_that("SPMD timeout class marks LL/LP CV bw opcodes as cv-regression", {
  timeout.class <- getFromNamespace(".npRmpi_spmd_timeout_class_from_opcode", "npRmpi")
  expect_identical(timeout.class("autodispatch.npregbw.cv_lllp"), "cv-regression")
  expect_identical(timeout.class("autodispatch.npscoefbw.cv_lllp"), "cv-regression")
  expect_identical(timeout.class("autodispatch.npplregbw.cv_lllp"), "cv-regression")
  expect_identical(timeout.class("autodispatch.npindexbw.core"), "cv-regression")
  expect_identical(timeout.class("autodispatch.npreg.core"), "default")
  expect_identical(timeout.class("autodispatch.npscoef.core"), "default")
  expect_identical(timeout.class("autodispatch.npplreg.core"), "default")
  expect_identical(timeout.class("autodispatch.npindex.core"), "default")
  expect_identical(timeout.class("autodispatch.npudens.core"), "default")
  expect_identical(timeout.class("autodispatch.npudist.core"), "default")
  expect_identical(timeout.class("autodispatch.npcdens.core"), "default")
  expect_identical(timeout.class("autodispatch.npcdist.core"), "default")
  expect_identical(timeout.class("autodispatch.npksum.core"), "default")
  expect_identical(timeout.class("autodispatch.npregiv.core"), "default")
  expect_identical(timeout.class("autodispatch.npcmstest.core"), "default")
  expect_identical(timeout.class("autodispatch.npindexbw"), "default")
})

test_that("SPMD npindexbw core opcode falls back when bandwidth.compute is FALSE", {
  opcode.fun <- getFromNamespace(".npRmpi_spmd_opcode_from_call", "npRmpi")
  mc.si <- quote(npindexbw(xdat = x, ydat = y, bandwidth.compute = FALSE))
  op.si <- opcode.fun(mc = mc.si, caller_env = environment())
  expect_identical(op.si, "autodispatch.npindexbw")
})

test_that("SPMD opcode selection tags density/distribution bw CV routes", {
  opcode.fun <- getFromNamespace(".npRmpi_spmd_opcode_from_call", "npRmpi")
  bw.flag <- TRUE

  mc.ud <- quote(npudensbw(dat = d, bandwidth.compute = TRUE))
  op.ud <- opcode.fun(mc = mc.ud, caller_env = environment())
  expect_identical(op.ud, "autodispatch.npudensbw.cv")

  mc.ui <- quote(npudistbw(dat = d))
  op.ui <- opcode.fun(mc = mc.ui, caller_env = environment())
  expect_identical(op.ui, "autodispatch.npudistbw.cv")

  mc.cd <- quote(npcdensbw(xdat = x, ydat = y, bandwidth.compute = bw.flag))
  op.cd <- opcode.fun(mc = mc.cd, caller_env = environment())
  expect_identical(op.cd, "autodispatch.npcdensbw.cv")

  mc.ci <- quote(npcdistbw(xdat = x, ydat = y, bandwidth.compute = TRUE))
  op.ci <- opcode.fun(mc = mc.ci, caller_env = environment())
  expect_identical(op.ci, "autodispatch.npcdistbw.cv")
})

test_that("SPMD density/distribution bw opcode falls back when bandwidth.compute is FALSE", {
  opcode.fun <- getFromNamespace(".npRmpi_spmd_opcode_from_call", "npRmpi")
  bw.flag <- FALSE

  mc.ud <- quote(npudensbw(dat = d, bandwidth.compute = FALSE))
  op.ud <- opcode.fun(mc = mc.ud, caller_env = environment())
  expect_identical(op.ud, "autodispatch.npudensbw")

  mc.ui <- quote(npudistbw(dat = d, bandwidth.compute = bw.flag))
  op.ui <- opcode.fun(mc = mc.ui, caller_env = environment())
  expect_identical(op.ui, "autodispatch.npudistbw")

  mc.cd <- quote(npcdensbw(xdat = x, ydat = y, bandwidth.compute = FALSE))
  op.cd <- opcode.fun(mc = mc.cd, caller_env = environment())
  expect_identical(op.cd, "autodispatch.npcdensbw")

  mc.ci <- quote(npcdistbw(xdat = x, ydat = y, bandwidth.compute = FALSE))
  op.ci <- opcode.fun(mc = mc.ci, caller_env = environment())
  expect_identical(op.ci, "autodispatch.npcdistbw")
})

test_that("SPMD timeout class marks density/distribution bw CV opcodes as cv-density", {
  timeout.class <- getFromNamespace(".npRmpi_spmd_timeout_class_from_opcode", "npRmpi")
  expect_identical(timeout.class("autodispatch.npudensbw.cv"), "cv-density")
  expect_identical(timeout.class("autodispatch.npudistbw.cv"), "cv-density")
  expect_identical(timeout.class("autodispatch.npcdensbw.cv"), "cv-density")
  expect_identical(timeout.class("autodispatch.npcdistbw.cv"), "cv-density")
})

test_that("SPMD locked core opcodes reject mismatched call heads", {
  make.env <- getFromNamespace(".npRmpi_spmd_make_envelope", "npRmpi")
  try.exec <- getFromNamespace(".npRmpi_spmd_try_execute_local", "npRmpi")
  set.seq <- getFromNamespace(".npRmpi_spmd_seq_set", "npRmpi")

  old.seq <- getOption("npRmpi.spmd.seq_id", 0L)
  on.exit(options(npRmpi.spmd.seq_id = old.seq), add = TRUE)
  set.seq(0L)

  env.reg <- make.env(opcode = "autodispatch.npregbw.cv_lllp", timeout_class = "cv-regression")
  bad.reg <- try.exec(env.reg, payload = list(call = quote(npudensbw(dat = d))), where = "unit locked opcode")
  expect_true(is.list(bad.reg) && !isTRUE(bad.reg$ok))
  expect_match(bad.reg$error, "restricted to")

  env.den <- make.env(opcode = "autodispatch.npudensbw.cv", timeout_class = "cv-density")
  bad.den <- try.exec(env.den, payload = list(call = quote(npregbw(y ~ x, regtype = "ll", bwmethod = "cv.ls"))),
                      where = "unit locked opcode")
  expect_true(is.list(bad.den) && !isTRUE(bad.den$ok))
  expect_match(bad.den$error, "restricted to")

  env.si <- make.env(opcode = "autodispatch.npindexbw.core", timeout_class = "cv-regression")
  bad.si <- try.exec(env.si, payload = list(call = quote(npplregbw(xdat = x, ydat = y, zdat = z))),
                     where = "unit locked opcode")
  expect_true(is.list(bad.si) && !isTRUE(bad.si$ok))
  expect_match(bad.si$error, "restricted to")

  env.si.est <- make.env(opcode = "autodispatch.npindex.core", timeout_class = "default")
  bad.si.est <- try.exec(env.si.est, payload = list(call = quote(npscoef(bws = bw))),
                         where = "unit locked opcode")
  expect_true(is.list(bad.si.est) && !isTRUE(bad.si.est$ok))
  expect_match(bad.si.est$error, "restricted to")

  env.reg.est <- make.env(opcode = "autodispatch.npreg.core", timeout_class = "default")
  bad.reg.est <- try.exec(env.reg.est, payload = list(call = quote(npudens(tdat = x, bws = bw))),
                          where = "unit locked opcode")
  expect_true(is.list(bad.reg.est) && !isTRUE(bad.reg.est$ok))
  expect_match(bad.reg.est$error, "restricted to")

  env.den.est <- make.env(opcode = "autodispatch.npcdens.core", timeout_class = "default")
  bad.den.est <- try.exec(env.den.est, payload = list(call = quote(npreg(bws = bw, gradients = FALSE))),
                          where = "unit locked opcode")
  expect_true(is.list(bad.den.est) && !isTRUE(bad.den.est$ok))
  expect_match(bad.den.est$error, "restricted to")

  env.qreg <- make.env(opcode = "autodispatch.npqreg.core", timeout_class = "default")
  bad.qreg <- try.exec(env.qreg, payload = list(call = quote(npconmode(bws = bw))),
                       where = "unit locked opcode")
  expect_true(is.list(bad.qreg) && !isTRUE(bad.qreg$ok))
  expect_match(bad.qreg$error, "restricted to")

  env.iv <- make.env(opcode = "autodispatch.npregiv.core", timeout_class = "default")
  bad.iv <- try.exec(env.iv, payload = list(call = quote(npregivderiv(y = y, z = z, w = w))),
                     where = "unit locked opcode")
  expect_true(is.list(bad.iv) && !isTRUE(bad.iv$ok))
  expect_match(bad.iv$error, "restricted to")
})

test_that("SPMD dynamic autodispatch opcode executes payload with ACK", {
  make.env <- getFromNamespace(".npRmpi_spmd_make_envelope", "npRmpi")
  try.exec <- getFromNamespace(".npRmpi_spmd_try_execute_local", "npRmpi")
  set.seq <- getFromNamespace(".npRmpi_spmd_seq_set", "npRmpi")

  old.seq <- getOption("npRmpi.spmd.seq_id", 0L)
  on.exit(options(npRmpi.spmd.seq_id = old.seq), add = TRUE)
  set.seq(0L)

  env <- make.env(opcode = "autodispatch.unit_test", timeout_class = "unit")
  payload <- list(call = quote(1 + 1))
  ans <- try.exec(envelope = env, payload = payload, where = "unit test")

  expect_true(is.list(ans) && isTRUE(ans$ok))
  expect_identical(ans$ack$status, "ACK")
  expect_identical(ans$ack$opcode, "autodispatch.unit_test")
  expect_identical(ans$result, 2)
})
