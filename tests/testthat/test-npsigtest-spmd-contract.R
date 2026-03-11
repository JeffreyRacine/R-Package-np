with_mpi_pool <- function(code) {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  force(code)
}

test_that("npsigtest maps to locked SPMD opcode", {
  opcode.fun <- getFromNamespace(".npRmpi_spmd_opcode_from_call", "npRmpi")

  mc.base <- quote(npsigtest(bws = bw, boot.num = 9))
  op.base <- opcode.fun(mc = mc.base, caller_env = environment())
  expect_identical(op.base, "autodispatch.npsigtest.core")

  mc.form <- quote(npsigtest.formula(y ~ x1 + x2, data = d, boot.num = 9))
  op.form <- opcode.fun(mc = mc.form, caller_env = environment())
  expect_identical(op.form, "autodispatch.npsigtest.core")
})

test_that("npsigtest locked opcode rejects mismatched call heads", {
  make.env <- getFromNamespace(".npRmpi_spmd_make_envelope", "npRmpi")
  try.exec <- getFromNamespace(".npRmpi_spmd_try_execute_local", "npRmpi")
  set.seq <- getFromNamespace(".npRmpi_spmd_seq_set", "npRmpi")

  with_mpi_pool({
    old.seq <- getOption("npRmpi.spmd.seq_id", 0L)
    on.exit(options(npRmpi.spmd.seq_id = old.seq), add = TRUE)
    set.seq(0L)

    env <- make.env(opcode = "autodispatch.npsigtest.core", timeout_class = "default")
    bad <- try.exec(env, payload = list(call = quote(npsymtest(data = y))), where = "unit locked opcode")
    expect_true(is.list(bad) && !isTRUE(bad$ok))
    expect_match(bad$error, "restricted to")
  })
})

test_that("npsigtest entrypoint orchestrates locally without whole-call autodispatch", {
  fn <- getFromNamespace("npsigtest", "npRmpi")
  body.txt <- paste(deparse(body(fn), width.cutoff = 500L), collapse = " ")
  expect_false(grepl("\\.npRmpi_autodispatch_call\\(", body.txt))
  expect_false(grepl("\\.npRmpi_manual_distributed_call\\(", body.txt))
})

test_that("npsigtest bandwidth extractor materializes xdat/ydat for formula bws", {
  ext <- getFromNamespace(".npRmpi_npsig_extract_xy_from_bws", "npRmpi")
  with_mpi_pool({
    set.seed(7)
    n <- 20L
    d <- data.frame(
      y = rnorm(n),
      x1 = runif(n),
      x2 = runif(n)
    )
    bw <- npregbw(y ~ x1 + x2, data = d, bws = c(0.2, 0.4), bandwidth.compute = FALSE)

    xy <- ext(bw)
    expect_true(is.list(xy))
    expect_true(is.data.frame(xy$xdat))
    expect_true(is.numeric(xy$ydat))
    expect_identical(nrow(xy$xdat), n)
    expect_identical(length(xy$ydat), n)
  })
})
