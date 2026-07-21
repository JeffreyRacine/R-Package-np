iv_summary_namespace <- function() {
  if ("npRmpi" %in% loadedNamespaces()) "npRmpi" else "np"
}

iv_summary_formatter <- function() {
  getFromNamespace(".np_iv_regression_summary_text", iv_summary_namespace())
}

test_that("IV summaries use canonical regression-estimator labels", {
  format <- iv_summary_formatter()
  lc <- format(list(effective = list(regtype = "lc", degree = 0L)),
               family = "npregiv")
  ll <- format(list(effective = list(regtype = "ll", degree = 1L)),
               family = "npregiv")
  lp <- format(list(effective = list(regtype = "lp", degree = 2L)),
               family = "npregiv")
  deriv.lp <- format(list(effective = list(regtype = "lp", degree = 2L)),
                     family = "npregivderiv")

  expect_identical(lc, "\nKernel Regression Estimator: Local-Constant")
  expect_identical(ll, "\nKernel Regression Estimator: Local-Linear")
  expect_identical(
    lp,
    paste0(
      "\nKernel Regression Estimator: ",
      "Local-Polynomial (Generalized basis; degree = 2)",
      "\nLP Basis Family: Generalized",
      "\nLP Basis Representation: ",
      "Shifted Legendre (degree-graded orthonormal)"
    )
  )
  expect_identical(
    deriv.lp,
    paste0(
      "\nKernel Regression Estimator: ",
      "Local-Polynomial (Generalized basis; degree = 2)",
      "\nLP Basis Family: Generalized",
      "\nLP Basis Representation: Raw"
    )
  )
})

test_that("IV smoothing text delegates to the canonical formatter", {
  ns <- iv_summary_namespace()
  format <- iv_summary_formatter()
  canonical <- getFromNamespace("genRegEstStr", ns)
  cases <- list(
    list(regtype = "lc", pregtype = "Local-Constant", degree = 0L),
    list(regtype = "ll", pregtype = "Local-Linear", degree = 1L),
    list(regtype = "lp", pregtype = "Local-Polynomial", degree = 2L)
  )

  for (case in cases) {
    expected <- canonical(c(case, list(
      basis = "glp", bernstein.basis = TRUE
    )))
    observed <- format(
      list(effective = case[c("regtype", "degree")]),
      family = "npregiv"
    )
    expect_identical(observed, expected)
  }
})

test_that("ordinary legacy p has a conservative canonical fallback", {
  format <- iv_summary_formatter()

  expect_identical(
    format(NULL, p = 0L, family = "npregiv"),
    "\nKernel Regression Estimator: Local-Constant"
  )
  expect_identical(
    format(NULL, p = 1L, family = "npregiv"),
    "\nKernel Regression Estimator: Local-Linear"
  )
  expect_match(
    format(NULL, p = 2L, family = "npregiv"),
    "Local-Polynomial \\(Generalized basis; degree = 2\\)"
  )
  expect_identical(format(NULL, p = 1L, family = "npregivderiv"), "")
})

test_that("malformed smoothing metadata fails closed", {
  format <- iv_summary_formatter()
  invalid <- list(
    list(effective = list(regtype = "ll", degree = 2L)),
    list(effective = list(regtype = "lc", degree = 1L)),
    list(effective = list(regtype = "bad", degree = 1L)),
    list(effective = list()),
    list(effective = "ll"),
    list(effective = list(regtype = "lp", degree = Inf)),
    list(effective = list(regtype = "lp", degree = c(1L, 2L))),
    list(effective = list(regtype = "lp", degree =
                           as.double(.Machine$integer.max) + 1))
  )
  expect_true(all(vapply(
    invalid,
    function(x) identical(format(x, p = 1L, family = "npregiv"), ""),
    logical(1)
  )))
  expect_identical(format(NULL, p = NA_real_, family = "npregiv"), "")
  expect_identical(format(NULL, p = 1.5, family = "npregiv"), "")
  expect_error(format(NULL, family = "unknown"), "arg")
})

test_that("LP representation labels remain tied to their engines", {
  ns <- iv_summary_namespace()
  w.lp <- getFromNamespace("W.lp", ns)
  rbandwidth <- getFromNamespace("rbandwidth", ns)

  expect_identical(eval(formals(w.lp)$basis),
                   c("glp", "additive", "tensor"))
  expect_identical(eval(formals(w.lp)$bernstein.basis), TRUE)
  expect_identical(eval(formals(rbandwidth)$basis),
                   c("glp", "additive", "tensor"))
  expect_identical(eval(formals(rbandwidth)$bernstein.basis), FALSE)
})

test_that("public IV summaries replace only the smoothing block", {
  if (identical(iv_summary_namespace(), "npRmpi"))
    skip("installed MPI summary sentinels own public-fit coverage")

  set.seed(2718)
  n <- 22L
  w <- rnorm(n)
  v <- rnorm(n, sd = 0.18)
  z <- 0.55 * w + v
  y <- sin(z) - 0.3 * v + rnorm(n, sd = 0.04)
  args <- list(y = y, z = z, w = w, nmulti = 1L, iterate.max = 2L,
               stop.on.increase = FALSE, random.seed = 2718L)
  ordinary <- suppressWarnings(suppressMessages(do.call(
    npregiv, c(args, list(regtype = "ll"))
  )))
  derivative <- suppressWarnings(suppressMessages(do.call(
    npregivderiv,
    c(args, list(iterate.break = FALSE, regtype = "ll"))
  )))

  for (fit in list(ordinary, derivative)) {
    object <- summary(fit)
    output <- capture.output(visible <- withVisible(print(object)))
    expect_false(visible$visible)
    expect_identical(visible$value, object)
    expect_identical(
      grep("^Kernel Regression Estimator:", output, value = TRUE),
      "Kernel Regression Estimator: Local-Linear"
    )
    expect_false(any(grepl("^Local polynomial order \\(p\\):", output)))
    expect_false(any(grepl("^Regression smoothing type:", output)))
    expect_false(any(grepl("^Regression smoothing degree:", output)))
    expect_true(any(grepl("^Training observations:", output)))
    expect_true(any(grepl("^Regularization method:", output)))
  }

  malformed <- summary(derivative)
  malformed$smoothing.spec <- list(
    effective = list(regtype = "ll", degree = 2L)
  )
  malformed.output <- capture.output(print(malformed))
  expect_false(any(grepl("^Kernel Regression Estimator:", malformed.output)))
  expect_true(any(grepl("^Training observations:", malformed.output)))
  expect_true(any(grepl("^Regularization method:", malformed.output)))
})
