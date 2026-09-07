# Single-index inference uses the regression engine on the fitted scalar index.
# These helpers do not select bandwidths or resample beta.
.np_index_kernel_sum <- function(..., bwtype) {
  if (is.null(.np_progress_runtime$fit_forward))
    return(npksum(..., bwtype = bwtype,
                  bandwidth.divide = identical(bwtype, "adaptive_nn")))
  # The parent reports activity, not a percentage across heterogeneous calls.
  # Enable the existing native row callbacks only during this scoped activity.
  args <- list(...)
  total <- max(NROW(args[["txdat"]]), NROW(args[["exdat"]]))
  .np_with_compiled_fit_progress(
    label = "Fitting single-index model", total = total,
    expr = npksum(..., bwtype = bwtype,
                  bandwidth.divide = identical(bwtype, "adaptive_nn")))
}

.np_index_asymptotic_outputs <- function(fit, beta, gradients = FALSE) {
  out <- list(merr = as.double(fit$merr))
  if (gradients) {
    out$gerr <- as.vector(fit$gerr[, 1L]) %o% abs(as.vector(beta))
  }
  out
}

.np_index_covariance_inverse <- function(information) {
  factor <- tryCatch(chol(information), error = function(e) {
    stop(paste0("npindex(): asymptotic coefficient covariance could not be computed: ",
                "the free-coefficient information matrix is singular or not positive definite. ",
                "Original factorization error: ", conditionMessage(e)),
         call. = FALSE)
  })
  chol2inv(factor)
}

.np_index_refit_hint <- function(expr, gradients = FALSE, se = FALSE) {
  object <- if (is.symbol(expr) && nzchar(as.character(expr))) {
    paste(deparse(expr, backtick = TRUE), collapse = "")
  } else {
    "object"
  }
  paste0("Refit without repeating bandwidth search: npindex(bws = ",
         object, "$bws, gradients = ", if (gradients) "TRUE" else "FALSE",
         ", se = ", if (se) "TRUE" else "FALSE", ").",
         if (identical(object, "object"))
           " Replace 'object' with your fitted model." else "")
}
