test_that("transport binds the already selected method without redispatch", {
  ns <- asNamespace("npRmpi")
  scope <- new.env(parent = ns)
  execute <- get(".npRmpi_spmd_eval_payload", ns)
  environment(execute) <- scope
  scope$.npRmpi_lease_prepare_local <- function(value, publication) list(kind = "none")
  scope$.npRmpi_eval_scmd <- function(call, envir) {
    expect_identical(call[[1L]], as.name("npudens"))
    expect_identical(get("npudens", envir, inherits = FALSE),
                     get("npudens.bandwidth", ns))
    37L
  }
  payload <- list(call = quote(npudens(bws = x)),
                  call.owner = "npudens.bandwidth", call.binding = "npudens",
                  publication = list(kind = "none"))
  expect_identical(execute(payload, NULL)$result, 37L)
  expect_false(isTRUE(getOption("npRmpi.autodispatch.context", FALSE)))
  expect_false(isTRUE(getOption("npRmpi.autodispatch.disable", FALSE)))
  payload$call.binding <- "some_other_name"
  expect_error(execute(payload, NULL), "invalid materialized package call binding")
  direct <- get(".npscoefbw_nomad_context_prepare", ns)
  payload$call <- as.call(list(direct))
  payload$call.owner <- payload$call.binding <- ".npscoefbw_nomad_context_prepare"
  scope$.npRmpi_eval_scmd <- function(call, envir) {
    expect_identical(call[[1L]], direct)
    41L
  }
  expect_identical(execute(payload, NULL)$result, 41L)
  payload$call[[1L]] <- function() NULL
  expect_error(execute(payload, NULL), "invalid materialized package call binding")
})
