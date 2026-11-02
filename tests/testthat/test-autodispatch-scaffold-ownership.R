test_that("package command temporaries do not occupy worker global bindings", {
  skip_on_cran()
  probe <- quote(local({
    library(npRmpi)
    npRmpi.init(nslaves = 1L, quiet = TRUE)
    on.exit(npRmpi.quit(force = TRUE), add = TRUE)
    options(np.messages = FALSE)
    ns <- asNamespace("npRmpi")
    rexec <- get("mpi.remote.exec", ns)
    for (collision in c(TRUE, FALSE)) {
      if (collision) {
        rexec(local({
          for (key in c("ans", "exec.step", "run", "fn",
                           "result", "value", "old.ctx", "old.disable"))
            assign(key, paste0("USER_", key), .GlobalEnv)
        }), ret = FALSE)
      } else {
        rexec(local({
          keys <- intersect(c("ans", "exec.step", "run", "fn",
                           "result", "value", "old.ctx", "old.disable"),
                            ls(.GlobalEnv, all.names = TRUE))
          if (length(keys)) rm(list = keys, envir = .GlobalEnv)
        }), ret = FALSE)
      }
      get(".npRmpi_session_reset_spmd_state_allranks", ns)()
      get(".npRmpi_spmd_tiny_smoke", ns)()
      d <- data.frame(x = seq(-1, 1, length.out = 31L))
      d$y <- d$x + sin(seq_len(31L)) * 0.2
      b <- npregbw(y ~ x, data = d, bws = 0.3, bandwidth.compute = FALSE)
      fit <- npreg(bws = b)
      stopifnot(all(is.finite(fitted(fit))))
      H <- npreghat(b, txdat = d["x"])
      stopifnot(identical(dim(H), c(31L, 31L)), all(is.finite(H)))
      failed <- try(npreghat(b, txdat = d["x"], output = "apply"), silent = TRUE)
      stopifnot(inherits(failed, "try-error"))
      stopifnot(identical(H, npreghat(b, txdat = d["x"])))
      rm(b, fit)
      invisible(gc())
      get(".npRmpi_lease_drain_pending", ns)()
      seen <- rexec(local({
        vapply(c("ans", "exec.step", "run", "fn",
                           "result", "value", "old.ctx", "old.disable"), function(key) {
          if (!exists(key, .GlobalEnv, inherits = FALSE)) return("absent")
          value <- get(key, .GlobalEnv)
          if (identical(value, paste0("USER_", key))) "sentinel" else "overwritten"
        }, character(1))
      }), simplify = FALSE)
      stopifnot(all(unlist(seen) == if (collision) "sentinel" else "absent"))
    }
    mpi.bcast.cmd(user_broadcast_probe <- 73L, caller.execute = FALSE)
    stopifnot(all(unlist(rexec(get("user_broadcast_probe", .GlobalEnv),
                               simplify = FALSE)) == 73L))
    cat("SCAFFOLD_OWNERSHIP_OK\n")
  }))
  result <- npRmpi_run_isolated_contract(deparse(probe, width.cutoff = 500L),
                                         marker = "SCAFFOLD_OWNERSHIP_OK", timeout = 45L)
  if (is.null(result)) skip("installed npRmpi unavailable for isolated contract")
  expect_identical(result$status, 0L, info = paste(result$output, collapse = "\n"))
  expect_true(result$witnessed)
})
