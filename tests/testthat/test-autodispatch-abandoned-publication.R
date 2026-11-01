test_that("acknowledged publications are retired when no result is returned", {
  skip_on_cran()
  probe <- quote(local({
  library(npRmpi)
  npRmpi.init(nslaves = 1L, quiet = TRUE)
  on.exit(npRmpi.quit(force = TRUE), add = TRUE)
  options(np.messages = FALSE)
  ns <- asNamespace("npRmpi")
  replace <- function(name, value) {
    unlockBinding(name, ns)
    assign(name, value, ns)
    lockBinding(name, ns)
  }
  commit <- get(".npRmpi_lease_commit_local", ns)
  tag <- get(".npRmpi_autodispatch_tag_result", ns)
  on.exit(replace(".npRmpi_lease_commit_local", commit), add = TRUE)
  on.exit(replace(".npRmpi_autodispatch_tag_result", tag), add = TRUE)
  inventory <- function() {
    c(list(ls(get(".npRmpi_lease_state", ns)$registry, all.names = TRUE)),
      get("mpi.remote.exec", ns)(
        ls(get(".npRmpi_lease_state", asNamespace("npRmpi"))$registry,
           all.names = TRUE), simplify = FALSE))
  }
  x <- seq(-1, 1, length.out = 31L)
  y <- x + sin(seq_along(x)) * 0.2
  for (fault in c("interrupt", "tag")) {
    before <- inventory()
    if (fault == "interrupt") {
      replace(".npRmpi_lease_commit_local", function(plan) {
        ans <- commit(plan)
        get(".npRmpi_mpi_interrupt_set_pending", ns)()
        ans
      })
    } else {
      replace(".npRmpi_autodispatch_tag_result", function(...) {
        stop(structure(list(message = "tag failure", call = NULL, detail = 27L),
                       class = c("test_tag_error", "error", "condition")))
      })
    }
    err <- tryCatch(npregbw(xdat = x, ydat = y, bws = 0.3,
                           bandwidth.compute = FALSE),
                    error = identity, interrupt = identity)
    replace(".npRmpi_lease_commit_local", commit)
    replace(".npRmpi_autodispatch_tag_result", tag)
    stopifnot(inherits(err, if (fault == "interrupt") "interrupt" else "test_tag_error"))
    if (fault == "tag") stopifnot(identical(err$detail, 27L))
    stopifnot(identical(inventory(), before))
    stopifnot(all(is.finite(npudens(tdat = x, bws = 0.3)$dens)))
  }
  cat("ABANDONED_PUBLICATION_OK\n")
  }))
  result <- npRmpi_run_isolated_contract(deparse(probe, width.cutoff = 500L),
                                         marker = "ABANDONED_PUBLICATION_OK", timeout = 45L)
  if (is.null(result)) skip("installed npRmpi unavailable for isolated contract")
  expect_identical(result$status, 0L, info = paste(result$output, collapse = "\n"))
  expect_true(result$witnessed)
})
