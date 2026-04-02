run_wrapper_cmd_subprocess <- function(cmd, args = character(), timeout = 60L, env = character()) {
  out <- suppressWarnings(system2(cmd,
                                  args,
                                  stdout = TRUE,
                                  stderr = TRUE,
                                  timeout = timeout,
                                  env = env))
  status <- attr(out, "status")
  if (is.null(status))
    status <- 0L
  list(status = as.integer(status), output = out)
}

.wrapper_mpi_init_env_failure <- function(output) {
  any(grepl("OFI call ep_enable failed", output, fixed = TRUE)) ||
    any(grepl("Fatal error in internal_Init", output, fixed = TRUE)) ||
    any(grepl("MPI_Init", output, fixed = TRUE) & grepl("failed", output, ignore.case = TRUE))
}

wrapper_subprocess_env <- function(extra = character()) {
  npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1", extra))
}

wrapper_subprocess_lib <- function(env) {
  libs <- env[grepl("^R_LIBS=", env)]
  if (!length(libs))
    return(NULL)
  strsplit(sub("^R_LIBS=", "", libs[[1L]]), .Platform$path.sep, fixed = TRUE)[[1L]][1L]
}

wrapper_session_lines <- function(lines) {
  c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "ns <- asNamespace('npRmpi')",
    "mpi_send <- get('mpi.send', envir = ns, inherits = FALSE)",
    "mpi_recv <- get('mpi.recv', envir = ns, inherits = FALSE)",
    "mpi_sendrecv <- get('mpi.sendrecv', envir = ns, inherits = FALSE)",
    "mpi_sendrecv_replace <- get('mpi.sendrecv.replace', envir = ns, inherits = FALSE)",
    "mpi_reduce <- get('mpi.reduce', envir = ns, inherits = FALSE)",
    "mpi_allreduce <- get('mpi.allreduce', envir = ns, inherits = FALSE)",
    "mpi_gather <- get('mpi.gather', envir = ns, inherits = FALSE)",
    "mpi_gatherv <- get('mpi.gatherv', envir = ns, inherits = FALSE)",
    "mpi_allgather <- get('mpi.allgather', envir = ns, inherits = FALSE)",
    "mpi_allgatherv <- get('mpi.allgatherv', envir = ns, inherits = FALSE)",
    "mpi_scatter <- get('mpi.scatter', envir = ns, inherits = FALSE)",
    "mpi_scatterv <- get('mpi.scatterv', envir = ns, inherits = FALSE)",
    "mpi_isend <- get('mpi.isend', envir = ns, inherits = FALSE)",
    "mpi_irecv <- get('mpi.irecv', envir = ns, inherits = FALSE)",
    "mpi_isend_Robj <- get('mpi.isend.Robj', envir = ns, inherits = FALSE)",
    "npRmpi.init(nslaves = 1L, quiet = TRUE)",
    "on.exit(try(npRmpi.quit(force = TRUE), silent = TRUE), add = TRUE)",
    lines
  )
}

test_that("wrapper preserved self and collective routes stay green in session subprocess", {
  env <- wrapper_subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for wrapper contract")

  res <- npRmpi_run_rscript_subprocess(
    lines = wrapper_session_lines(c(
      "rank <- mpi.comm.rank(0L)",
      "stopifnot(identical(as.numeric(mpi_sendrecv(3.25, sendtype = 2L, dest = rank, sendtag = 31L, recvdata = double(1), recvtype = 2L, source = rank, recvtag = 31L, comm = 0L, status = 0L)), 3.25))",
      "stopifnot(identical(as.numeric(mpi.bcast(3.25, type = 5L, rank = 0L, comm = 0L, buffunit = 100L)), 3.25))",
      "stopifnot(identical(as.numeric(mpi_reduce(3.25, type = 2L, op = 'sum', dest = 0L, comm = 0L)), 3.25))",
      "stopifnot(identical(as.numeric(mpi_allreduce(3.25, type = 2L, op = 'sum', comm = 0L)), 3.25))",
      "stopifnot(identical(as.integer(mpi_gather(7L, type = 1L, rdata = integer(1), root = 0L, comm = 0L)), 7L))",
      "stopifnot(identical(as.integer(mpi_gatherv(as.raw(7), type = 4L, rdata = raw(1), rcounts = 1L, root = 0L, comm = 0L)), 7L))",
      "stopifnot(identical(as.integer(mpi_allgather(7L, type = 1L, rdata = integer(1), comm = 0L)), 7L))",
      "stopifnot(identical(as.integer(mpi_allgatherv(as.raw(7), type = 4L, rdata = raw(1), rcounts = 1L, comm = 0L)), 7L))",
      "stopifnot(identical(as.integer(mpi_scatter(7L, type = 1L, rdata = integer(1), root = 0L, comm = 0L)), 7L))",
      "stopifnot(identical(as.integer(mpi_scatterv(as.raw(7), scounts = 1L, type = 4L, rdata = raw(1), root = 0L, comm = 0L)), 7L))",
      "stopifnot(identical(as.numeric(mpi_sendrecv_replace(3.25, type = 5L, dest = rank, sendtag = 41L, source = rank, recvtag = 41L, comm = 0L, status = 0L)), 3.25))",
      "cat('WRAPPER_SESSION_POSITIVE_OK\\n')"
    )),
    timeout = 45L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("WRAPPER_SESSION_POSITIVE_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("wrapper quarantined routes fail fast with explicit errors in session subprocess", {
  env <- wrapper_subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for wrapper contract")

  res <- npRmpi_run_rscript_subprocess(
    lines = wrapper_session_lines(c(
      "msg_send <- tryCatch({ mpi_send(as.raw(1), type = 5L, dest = 0L, tag = 1L, comm = 0L); '' }, error = conditionMessage)",
      "msg_recv <- tryCatch({ mpi_recv(raw(1), type = 5L, source = 0L, tag = 1L, comm = 0L, status = 0L); '' }, error = conditionMessage)",
      "msg_bcast <- tryCatch({ mpi.bcast(as.raw(1), type = 6L, rank = 0L, comm = 0L); '' }, error = conditionMessage)",
      "msg_reduce <- tryCatch({ .Call('mpi_reduce', as.raw(1), 4L, 1L, 0L, 0L, PACKAGE = 'npRmpi'); '' }, error = conditionMessage)",
      "msg_allreduce <- tryCatch({ .Call('mpi_allreduce', as.raw(1), 4L, 1L, 0L, PACKAGE = 'npRmpi'); '' }, error = conditionMessage)",
      "msg_sendrecv <- tryCatch({ mpi_sendrecv(3.25, sendtype = 5L, dest = 0L, sendtag = 1L, recvdata = double(1), recvtype = 2L, source = 0L, recvtag = 1L, comm = 0L, status = 0L); '' }, error = conditionMessage)",
      "msg_replace <- tryCatch({ mpi_sendrecv_replace('x', type = 3L, dest = 0L, sendtag = 1L, source = 0L, recvtag = 1L, comm = 0L, status = 0L); '' }, error = conditionMessage)",
      "msg_gather <- tryCatch({ mpi_gather(as.raw(1), type = 5L, rdata = raw(1), root = 0L, comm = 0L); '' }, error = conditionMessage)",
      "msg_gatherv <- tryCatch({ mpi_gatherv(as.raw(1), type = 5L, rdata = raw(1), rcounts = 1L, root = 0L, comm = 0L); '' }, error = conditionMessage)",
      "msg_allgather <- tryCatch({ mpi_allgather(as.raw(1), type = 5L, rdata = raw(1), comm = 0L); '' }, error = conditionMessage)",
      "msg_allgatherv <- tryCatch({ mpi_allgatherv(as.raw(1), type = 5L, rdata = raw(1), rcounts = 1L, comm = 0L); '' }, error = conditionMessage)",
      "msg_gather_chr <- tryCatch({ mpi_gather('abc', type = 3L, rdata = string(4), root = 0L, comm = 0L); '' }, error = conditionMessage)",
      "msg_gatherv_chr <- tryCatch({ mpi_gatherv('abc', type = 3L, rdata = string(4), rcounts = 3L, root = 0L, comm = 0L); '' }, error = conditionMessage)",
      "msg_allgather_chr <- tryCatch({ mpi_allgather('abc', type = 3L, rdata = string(4), comm = 0L); '' }, error = conditionMessage)",
      "msg_allgatherv_chr <- tryCatch({ mpi_allgatherv('abc', type = 3L, rdata = string(4), rcounts = 3L, comm = 0L); '' }, error = conditionMessage)",
      "msg_scatter <- tryCatch({ mpi_scatter(as.raw(1), type = 5L, rdata = raw(1), root = 0L, comm = 0L); '' }, error = conditionMessage)",
      "msg_scatterv <- tryCatch({ mpi_scatterv(as.raw(1), scounts = 1L, type = 5L, rdata = raw(1), root = 0L, comm = 0L); '' }, error = conditionMessage)",
      "msg_isend <- tryCatch({ mpi_isend(1L, type = 1L, dest = 0L, tag = 1L, comm = 0L, request = 0L); '' }, error = conditionMessage)",
      "msg_irecv <- tryCatch({ mpi_irecv(integer(1), type = 1L, source = 0L, tag = 1L, comm = 0L, request = 0L); '' }, error = conditionMessage)",
      "msg_isend_robj <- tryCatch({ mpi_isend_Robj(list(a = 1L), dest = 0L, tag = 1L, comm = 0L, request = 0L); '' }, error = conditionMessage)",
      "stopifnot(identical(msg_send, 'mpi_send: unsupported type 5'))",
      "stopifnot(identical(msg_recv, 'mpi_recv: unsupported type 5'))",
      "stopifnot(identical(msg_bcast, 'mpi_bcast: unsupported type 6'))",
      "stopifnot(identical(msg_reduce, 'mpi_reduce: only integer and double types are supported'))",
      "stopifnot(identical(msg_allreduce, 'mpi_allreduce: only integer and double types are supported'))",
      "stopifnot(identical(msg_sendrecv, 'mpi_sendrecv: unsupported send type 5'))",
      "stopifnot(identical(msg_replace, 'mpi_sendrecv_replace: character type is unsupported'))",
      "stopifnot(identical(msg_gather, 'mpi_gather: unsupported type code; only types 1-4 are supported'))",
      "stopifnot(identical(msg_gatherv, 'mpi_gatherv: unsupported type code; only types 1-4 are supported'))",
      "stopifnot(identical(msg_allgather, 'mpi_allgather: unsupported type code; only types 1-4 are supported'))",
      "stopifnot(identical(msg_allgatherv, 'mpi_allgatherv: unsupported type code; only types 1-4 are supported'))",
      "stopifnot(identical(msg_gather_chr, 'mpi_gather: character collectives are unsupported; use mpi.gather.Robj() or raw serialization'))",
      "stopifnot(identical(msg_gatherv_chr, 'mpi_gatherv: character collectives are unsupported; use mpi.gather.Robj() or raw serialization'))",
      "stopifnot(identical(msg_allgather_chr, 'mpi_allgather: character collectives are unsupported; use mpi.allgather.Robj() or raw serialization'))",
      "stopifnot(identical(msg_allgatherv_chr, 'mpi_allgatherv: character collectives are unsupported; use mpi.allgather.Robj() or raw serialization'))",
      "stopifnot(identical(msg_scatter, 'mpi_scatter: unsupported type code; only types 1-4 are supported'))",
      "stopifnot(identical(msg_scatterv, 'mpi_scatterv: unsupported type code; only types 1-4 are supported'))",
      "stopifnot(identical(msg_isend, 'mpi.isend is temporarily unsupported in npRmpi; use blocking mpi.send() or mpi.send.Robj() instead'))",
      "stopifnot(identical(msg_irecv, 'mpi.irecv is temporarily unsupported in npRmpi; use blocking mpi.recv() or mpi.recv.Robj() instead'))",
      "stopifnot(identical(msg_isend_robj, 'mpi.isend.Robj is temporarily unsupported in npRmpi; use blocking mpi.send.Robj() instead'))",
      "cat('WRAPPER_SESSION_NEGATIVE_OK\\n')"
    )),
    timeout = 45L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("WRAPPER_SESSION_NEGATIVE_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("supported string sendrecv, broadcast, and scatter routes stay green in session subprocess", {
  env <- wrapper_subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for wrapper contract")

  res <- npRmpi_run_rscript_subprocess(
    lines = wrapper_session_lines(c(
      "rank <- mpi.comm.rank(0L)",
      "stopifnot(identical(as.character(mpi_sendrecv('abc', sendtype = 3L, dest = rank, sendtag = 31L, recvdata = string(4), recvtype = 3L, source = rank, recvtag = 31L, comm = 0L, status = 0L)), 'abc'))",
      "stopifnot(identical(as.character(mpi.bcast('abc', type = 3L, rank = 0L, comm = 0L)), 'abc'))",
      "stopifnot(identical(as.character(mpi_scatter('abc', type = 3L, rdata = string(4), root = 0L, comm = 0L)), 'abc'))",
      "stopifnot(identical(as.character(mpi_scatterv('abc', scounts = 3L, type = 3L, rdata = string(4), root = 0L, comm = 0L)), 'abc'))",
      "cat('WRAPPER_SESSION_STRING_OK\\n')"
    )),
    timeout = 45L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("WRAPPER_SESSION_STRING_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("wrapper attach smoke stays green under mpiexec when enabled", {
  skip_if_not(identical(Sys.getenv("NP_RMPI_ENABLE_ATTACH_TEST"), "1"),
              "set NP_RMPI_ENABLE_ATTACH_TEST=1 to run attach-mode wrapper smoke")

  mpiexec <- Sys.which("mpiexec")
  skip_if(!nzchar(mpiexec), "mpiexec unavailable")

  env <- wrapper_subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for wrapper contract")

  script <- tempfile("npRmpi-wrapper-attach-", fileext = ".R")
  on.exit(unlink(script), add = TRUE)
  writeLines(c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "ns <- asNamespace('npRmpi')",
    "mpi_sendrecv_replace <- get('mpi.sendrecv.replace', envir = ns, inherits = FALSE)",
    "is.master <- isTRUE(npRmpi.init(mode = 'attach', quiet = TRUE, autodispatch = TRUE))",
    "on.exit({",
    "  try(npRmpi.quit(mode = 'attach'), silent = TRUE)",
    "  if (isTRUE(is.master)) try(Rmpi::mpi.quit(), silent = TRUE)",
    "}, add = TRUE)",
    "if (isTRUE(is.master)) {",
    "  stopifnot(identical(as.numeric(mpi_sendrecv_replace(3.25, type = 5L, dest = 0L, sendtag = 41L, source = 0L, recvtag = 41L, comm = 0L, status = 0L)), 3.25))",
    "  stopifnot(identical(as.numeric(mpi.bcast(3.25, type = 5L, rank = 0L, comm = 0L, buffunit = 100L)), 3.25))",
    "  cat('WRAPPER_ATTACH_OK\\n')",
    "}"
  ), script, useBytes = TRUE)

  run_once <- function(iface) {
    run_wrapper_cmd_subprocess(
      mpiexec,
      args = c("-n", "2", file.path(R.home("bin"), "Rscript"), "--no-save", script),
      timeout = 60L,
      env = c(
        env,
        "R_PROFILE_USER=",
        "R_PROFILE=",
        sprintf("FI_TCP_IFACE=%s", iface),
        "FI_PROVIDER=tcp",
        sprintf("FI_SOCKETS_IFACE=%s", iface)
      )
    )
  }

  res <- run_once("en0")
  if (res$status != 0L)
    res <- run_once("lo0")

  if (res$status != 0L && .wrapper_mpi_init_env_failure(res$output))
    skip("MPI runtime interface unavailable in this environment for attach-mode wrapper smoke")

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("WRAPPER_ATTACH_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("wrapper profile smoke stays green under mpiexec when enabled", {
  skip_if_not(identical(Sys.getenv("NP_RMPI_ENABLE_PROFILE_TEST"), "1"),
              "set NP_RMPI_ENABLE_PROFILE_TEST=1 to run profile-mode wrapper smoke")

  mpiexec <- Sys.which("mpiexec")
  skip_if(!nzchar(mpiexec), "mpiexec unavailable")

  env <- wrapper_subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for wrapper contract")

  lib.path <- wrapper_subprocess_lib(env)
  skip_if(is.null(lib.path), "local npRmpi install unavailable for wrapper contract")
  profile.path <- file.path(lib.path, "npRmpi", "Rprofile")
  skip_if(!file.exists(profile.path), "npRmpi profile template unavailable in subprocess lib")

  script <- tempfile("npRmpi-wrapper-profile-", fileext = ".R")
  on.exit(unlink(script), add = TRUE)
  writeLines(c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "if (mpi.comm.rank(0L) == 0L) {",
    "  ns <- asNamespace('npRmpi')",
    "  mpi_sendrecv_replace <- get('mpi.sendrecv.replace', envir = ns, inherits = FALSE)",
    "  mpi.bcast.cmd(np.mpi.initialize(), caller.execute = TRUE)",
    "  stopifnot(identical(as.numeric(mpi_sendrecv_replace(3.25, type = 5L, dest = 0L, sendtag = 41L, source = 0L, recvtag = 41L, comm = 0L, status = 0L)), 3.25))",
    "  stopifnot(identical(as.numeric(mpi.bcast(3.25, type = 5L, rank = 0L, comm = 0L, buffunit = 100L)), 3.25))",
    "  cat('WRAPPER_PROFILE_OK\\n')",
    "  mpi.bcast.cmd(mpi.quit(), caller.execute = TRUE)",
    "}"
  ), script, useBytes = TRUE)

  run_once <- function(iface) {
    run_wrapper_cmd_subprocess(
      mpiexec,
      args = c("-n", "2", file.path(R.home("bin"), "Rscript"), "--no-save", script),
      timeout = 60L,
      env = c(
        env,
        sprintf("R_PROFILE_USER=%s", profile.path),
        "R_PROFILE=",
        "NP_RMPI_PROFILE_RECV_TIMEOUT_SEC=60",
        sprintf("FI_TCP_IFACE=%s", iface),
        "FI_PROVIDER=tcp",
        sprintf("FI_SOCKETS_IFACE=%s", iface)
      )
    )
  }

  res <- run_once("en0")
  if (res$status != 0L)
    res <- run_once("lo0")

  if (res$status != 0L && .wrapper_mpi_init_env_failure(res$output))
    skip("MPI runtime interface unavailable in this environment for profile-mode wrapper smoke")

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("WRAPPER_PROFILE_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("supported string send/recv and bcast stay green under profile mpiexec when enabled", {
  skip_if_not(identical(Sys.getenv("NP_RMPI_ENABLE_PROFILE_TEST"), "1"),
              "set NP_RMPI_ENABLE_PROFILE_TEST=1 to run profile-mode wrapper smoke")

  mpiexec <- Sys.which("mpiexec")
  skip_if(!nzchar(mpiexec), "mpiexec unavailable")

  env <- wrapper_subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for wrapper contract")

  lib.path <- wrapper_subprocess_lib(env)
  skip_if(is.null(lib.path), "local npRmpi install unavailable for wrapper contract")
  profile.path <- file.path(lib.path, "npRmpi", "Rprofile")
  skip_if(!file.exists(profile.path), "npRmpi profile template unavailable in subprocess lib")

  out.dir <- tempfile("npRmpi-wrapper-profile-string-")
  dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(out.dir, recursive = TRUE, force = TRUE), add = TRUE)
  recv.path <- file.path(out.dir, "recv_rank1.txt")
  bcast0.path <- file.path(out.dir, "bcast_rank0.txt")
  bcast1.path <- file.path(out.dir, "bcast_rank1.txt")

  script <- tempfile("npRmpi-wrapper-profile-string-", fileext = ".R")
  on.exit(unlink(script), add = TRUE)
  writeLines(c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "if (mpi.comm.rank(0L) == 0L) {",
    "  mpi.bcast.cmd(np.mpi.initialize(), caller.execute = TRUE)",
    "  mpi.bcast.cmd({",
    "    ns <- asNamespace('npRmpi')",
    "    mpi_send <- get('mpi.send', envir = ns, inherits = FALSE)",
    "    mpi_recv <- get('mpi.recv', envir = ns, inherits = FALSE)",
    "    rank <- mpi.comm.rank(0L)",
    "    if (rank == 0L) {",
    "      mpi_send('abc', type = 3L, dest = 1L, tag = 11L, comm = 0L)",
    "    } else {",
    "      out <- mpi_recv(string(4), type = 3L, source = 0L, tag = 11L, comm = 0L, status = 0L)",
    sprintf("      writeLines(out, %s)", deparse(recv.path)),
    "    }",
    "  }, caller.execute = TRUE)",
    "  mpi.bcast.cmd({",
    "    rank <- mpi.comm.rank(0L)",
    "    out <- mpi.bcast(if (rank == 0L) 'abc' else string(4), type = 3L, rank = 0L, comm = 0L)",
    "    if (rank == 0L) {",
    sprintf("      writeLines(out, %s)", deparse(bcast0.path)),
    "    } else {",
    sprintf("      writeLines(out, %s)", deparse(bcast1.path)),
    "    }",
    "  }, caller.execute = TRUE)",
    "  mpi.bcast.cmd(mpi.quit(), caller.execute = TRUE)",
    "}"
  ), script, useBytes = TRUE)

  run_once <- function(iface) {
    run_wrapper_cmd_subprocess(
      mpiexec,
      args = c("-n", "2", file.path(R.home("bin"), "Rscript"), "--no-save", script),
      timeout = 60L,
      env = c(
        env,
        sprintf("R_PROFILE_USER=%s", profile.path),
        "R_PROFILE=",
        "NP_RMPI_PROFILE_RECV_TIMEOUT_SEC=60",
        sprintf("FI_TCP_IFACE=%s", iface),
        "FI_PROVIDER=tcp",
        sprintf("FI_SOCKETS_IFACE=%s", iface)
      )
    )
  }

  res <- run_once("en0")
  if (res$status != 0L)
    res <- run_once("lo0")

  if (res$status != 0L && .wrapper_mpi_init_env_failure(res$output))
    skip("MPI runtime interface unavailable in this environment for profile-mode wrapper smoke")

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(file.exists(recv.path))
  expect_equal(readLines(recv.path, warn = FALSE), "abc")
  expect_equal(readLines(bcast0.path, warn = FALSE), "abc")
  expect_equal(readLines(bcast1.path, warn = FALSE), "abc")
})
