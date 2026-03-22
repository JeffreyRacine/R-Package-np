gsl_bs_nprmpi <- getFromNamespace("gsl.bs", "npRmpi")
predict_gsl_bs_nprmpi <- getFromNamespace("predict.gsl.bs", "npRmpi")

capture_gsl_bs_eval <- function(expr) {
  warnings <- character()

  value <- tryCatch(
    withCallingHandlers(
      expr,
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      structure(
        list(
          message = conditionMessage(e),
          class = class(e)
        ),
        class = "captured_error"
      )
    }
  )

  if (inherits(value, "captured_error")) {
    return(list(
      ok = FALSE,
      warnings = warnings,
      error_message = value$message,
      error_class = value$class
    ))
  }

  list(
    ok = TRUE,
    warnings = warnings,
    value = value
  )
}

predict_contract_cases <- list(
  uniform_intercept = list(
    x = c(-1, -0.5, 0, 0.5, 1),
    degree = 3,
    nbreak = 4,
    deriv = 0,
    intercept = TRUE,
    x.min = NULL,
    x.max = NULL,
    knots = NULL
  ),
  deriv_knots_no_intercept = list(
    x = c(-1, -0.5, 0, 0.5, 1),
    degree = 3,
    nbreak = 4,
    deriv = 1,
    intercept = FALSE,
    x.min = -1,
    x.max = 1,
    knots = c(-1, -0.3, 0.2, 1)
  )
)

predict_contract_newx_cases <- list(
  newx_null = NULL,
  numeric_vector = c(-1, -0.25, 0.25, 1),
  integer_vector = as.integer(c(-1, 0, 1)),
  scalar_numeric = 0.2,
  numeric_matrix = matrix(c(-1, -0.25, 0.25, 1), ncol = 1),
  integer_matrix = matrix(as.integer(c(-1, 0, 1)), ncol = 1),
  named_vector = setNames(c(-1, 0, 1), c("a", "b", "c")),
  duplicates_outside = c(-2, -1, -1, 0, 1, 2),
  zero_length = numeric(0),
  na_nan_inf = c(NA_real_, NaN, Inf, -Inf, 0)
)

test_that("npRmpi gsl.bs uses package-specific first class with compatibility superclass", {
  obj <- gsl_bs_nprmpi(c(-1, -0.5, 0, 0.5, 1), degree = 3, nbreak = 4, intercept = TRUE)

  expect_identical(class(obj), c("npRmpi_gsl.bs", "gsl.bs", "matrix"))
  expect_true(inherits(obj, "gsl.bs"))
  expect_true(inherits(obj, "matrix"))
  expect_true(is.function(getS3method("predict", "npRmpi_gsl.bs")))
})

test_that("npRmpi gsl.bs predict dispatch matches legacy wrapper across supported newx forms", {
  method <- getS3method("predict", "npRmpi_gsl.bs")

  for (case_name in names(predict_contract_cases)) {
    args_case <- predict_contract_cases[[case_name]]
    obj <- do.call(gsl_bs_nprmpi, args_case)

    expect_identical(class(obj), c("npRmpi_gsl.bs", "gsl.bs", "matrix"))

    for (newx_name in names(predict_contract_newx_cases)) {
      newx <- predict_contract_newx_cases[[newx_name]]

      dispatch_res <- capture_gsl_bs_eval(predict(obj, newx = newx))
      method_res <- capture_gsl_bs_eval(method(obj, newx = newx))
      legacy_res <- capture_gsl_bs_eval(predict_gsl_bs_nprmpi(obj, newx = newx))

      expect_identical(
        dispatch_res,
        method_res,
        info = sprintf("%s/%s dispatch vs method", case_name, newx_name)
      )
      expect_identical(
        dispatch_res,
        legacy_res,
        info = sprintf("%s/%s dispatch vs legacy", case_name, newx_name)
      )

      if (dispatch_res$ok) {
        if (is.null(newx)) {
          expected <- obj
        } else {
          expected <- gsl_bs_nprmpi(
            x = as.numeric(newx),
            degree = attr(obj, "degree"),
            nbreak = attr(obj, "nbreak"),
            deriv = attr(obj, "deriv"),
            x.min = attr(obj, "x.min"),
            x.max = attr(obj, "x.max"),
            intercept = attr(obj, "intercept"),
            knots = attr(obj, "knots")
          )
        }

        attr(expected, "newx") <- if (is.null(newx)) NULL else as.numeric(newx)
        attr(expected, "newx.trimmed") <- NULL

        expect_identical(
          dispatch_res$value,
          expected,
          info = sprintf("%s/%s semantic parity", case_name, newx_name)
        )
      }
    }
  }
})
