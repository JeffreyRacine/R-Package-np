capture_np_error <- function(expr) {
  tryCatch(force(expr), error = identity)
}

test_that("conditional beta errors name the failing public kernel side", {
  x <- data.frame(x = seq(0.05, 0.95, length.out = 24L))
  y <- data.frame(y = seq(0.08, 0.92, length.out = 24L))
  common <- list(
    xdat = x, ydat = y, bws = c(0.2, 0.2),
    bandwidth.compute = FALSE, regtype = "lp", degree = 1L
  )

  x.invalid <- c(common, list(
    cxkertype = "beta", cykertype = "beta", cykerbound = "range"
  ))
  y.invalid <- c(common, list(
    cxkertype = "gaussian", cykertype = "beta", cxkerbound = "range"
  ))

  expected.x <- paste(
    "conditional density explanatory beta kernel requires",
    "cxkerbound = \"fixed\" or \"range\" with finite cxkerlb and cxkerub"
  )
  expected.y <- paste(
    "conditional distribution dependent beta kernel requires",
    "cykerbound = \"fixed\" or \"range\" with finite cykerlb and cykerub"
  )

  x.condition <- capture_np_error(do.call(npcdensbw, x.invalid))
  y.condition <- capture_np_error(do.call(npcdistbw, y.invalid))

  expect_s3_class(x.condition, "simpleError")
  expect_s3_class(y.condition, "simpleError")
  expect_identical(conditionMessage(x.condition), expected.x)
  expect_identical(conditionMessage(y.condition), expected.y)
  expect_null(conditionCall(x.condition))
  expect_null(conditionCall(y.condition))

  formula.data <- data.frame(y = y$y, x = x$x)
  formula.x <- tryCatch(do.call(npcdens, list(
    formula = y ~ x, data = formula.data,
    cxkertype = "beta", cykertype = "beta", cykerbound = "range",
    nomad = TRUE, bwmethod = "cv.ls"
  )), error = identity)
  formula.y <- tryCatch(do.call(npcdist, list(
    formula = y ~ x, data = formula.data,
    cxkertype = "gaussian", cykertype = "beta", cxkerbound = "range",
    nomad = TRUE, bwmethod = "cv.ls"
  )), error = identity)
  expect_identical(conditionMessage(formula.x), expected.x)
  expect_identical(conditionMessage(formula.y), expected.y)
})

test_that("conditional beta validation retains X-before-Y failure precedence", {
  x <- data.frame(x = seq(0.05, 0.95, length.out = 24L))
  y <- data.frame(y = seq(0.08, 0.92, length.out = 24L))

  condition <- capture_np_error(npcdensbw(
    xdat = x, ydat = y, bws = c(0.2, 0.2),
    bandwidth.compute = FALSE, regtype = "lp", degree = 1L,
    cxkertype = "beta", cykertype = "beta"
  ))

  expect_identical(
    conditionMessage(condition),
    paste(
      "conditional density explanatory beta kernel requires",
      "cxkerbound = \"fixed\" or \"range\" with finite cxkerlb and cxkerub"
    )
  )
})

conditional_kernel_lines <- function(object, method = c("summary", "print")) {
  method <- match.arg(method)
  output <- capture.output(if (identical(method, "summary")) {
    summary(object)
  } else {
    print(object)
  })
  trimws(grep("^\\s*Continuous Kernel Type", output, value = TRUE))
}

make_conditional_kernel_bw <- function(kind = c("density", "distribution"),
                                       cxkertype = "beta",
                                       cykertype = "beta",
                                       cxkerorder = 4L,
                                       cykerorder = 6L,
                                       cxkerbound = "range",
                                       cykerbound = "range") {
  kind <- match.arg(kind)
  x <- data.frame(x = seq(0.05, 0.95, length.out = 24L))
  y <- data.frame(y = seq(0.08, 0.92, length.out = 24L))
  do.call(
    if (identical(kind, "density")) npcdensbw else npcdistbw,
    list(
      xdat = x, ydat = y, bws = c(0.2, 0.2),
      bandwidth.compute = FALSE, regtype = "lp", degree = 1L,
      cxkertype = cxkertype, cykertype = cykertype,
      cxkerorder = cxkerorder, cykerorder = cykerorder,
      cxkerbound = cxkerbound, cykerbound = cykerbound,
      cykerlb = if (identical(cykerbound, "fixed")) 0 else NULL,
      cykerub = if (identical(cykerbound, "fixed")) 1 else NULL
    )
  )
}

test_that("conditional summaries always identify X and Y continuous kernels", {
  expected.same <- c(
    "Continuous Kernel Type (Exp. Var.): Fourth-Order Beta associated (bounded/range)",
    "Continuous Kernel Type (Dep. Var.): Fourth-Order Beta associated (bounded/range)"
  )
  expected.different <- c(
    "Continuous Kernel Type (Exp. Var.): Fourth-Order Beta associated (bounded/range)",
    "Continuous Kernel Type (Dep. Var.): Sixth-Order Beta associated (bounded/range)"
  )

  for (kind in c("density", "distribution")) {
    same <- make_conditional_kernel_bw(
      kind, cxkerorder = 4L, cykerorder = 4L
    )
    different <- make_conditional_kernel_bw(kind)

    expect_identical(conditional_kernel_lines(same), expected.same)
    expect_identical(conditional_kernel_lines(same, "print"), expected.same)
    expect_identical(conditional_kernel_lines(different), expected.different)
    expect_identical(
      conditional_kernel_lines(different, "print"), expected.different
    )
  }
})

test_that("conditional summaries retain side-specific boundary descriptions", {
  object <- make_conditional_kernel_bw(
    "density", cxkerorder = 4L, cykerorder = 4L,
    cxkerbound = "range", cykerbound = "fixed"
  )

  expect_identical(
    conditional_kernel_lines(object),
    c(
      "Continuous Kernel Type (Exp. Var.): Fourth-Order Beta associated (bounded/range)",
      "Continuous Kernel Type (Dep. Var.): Fourth-Order Beta associated (bounded/fixed)"
    )
  )
})

test_that("conditional summaries print only active continuous sides", {
  x <- data.frame(x = seq(0.05, 0.95, length.out = 24L))
  y <- data.frame(y = factor(rep(c("a", "b"), 12L)))
  object <- npcdensbw(
    xdat = x, ydat = y, bws = c(0.2, 0.2), bandwidth.compute = FALSE,
    regtype = "lp", degree = 1L,
    cxkertype = "beta", cxkerbound = "range"
  )

  expected <- paste(
    "Continuous Kernel Type (Exp. Var.):",
    "Second-Order Beta associated (bounded/range)"
  )
  expect_identical(conditional_kernel_lines(object), expected)
  expect_identical(conditional_kernel_lines(object, "print"), expected)
})

test_that("inherited conditional estimators use the canonical kernel lines", {
  distribution.bw <- make_conditional_kernel_bw(
    "distribution", cykerorder = 2L
  )
  quantile <- npqreg(bws = distribution.bw, tau = 0.5)

  x <- data.frame(x = seq(0.05, 0.95, length.out = 24L))
  y <- data.frame(y = factor(rep(c("a", "b"), 12L)))
  mode.bw <- npcdensbw(
    xdat = x, ydat = y, bws = c(0.2, 0.2), bandwidth.compute = FALSE,
    regtype = "lp", degree = 1L,
    cxkertype = "beta", cxkerbound = "range"
  )
  mode <- npconmode(bws = mode.bw)

  expect_identical(
    conditional_kernel_lines(quantile),
    c(
      "Continuous Kernel Type (Exp. Var.): Fourth-Order Beta associated (bounded/range)",
      "Continuous Kernel Type (Dep. Var.): Second-Order Beta associated (bounded/range)"
    )
  )
  expect_identical(
    conditional_kernel_lines(mode),
    paste(
      "Continuous Kernel Type (Exp. Var.):",
      "Second-Order Beta associated (bounded/range)"
    )
  )
})

test_that("unconditional summaries retain the generic cker display", {
  x <- data.frame(x = seq(0.05, 0.95, length.out = 24L))
  object <- npudensbw(
    dat = x, bws = 0.2, bandwidth.compute = FALSE,
    ckertype = "beta", ckerorder = 4L, ckerbound = "range"
  )

  expect_identical(
    conditional_kernel_lines(object),
    "Continuous Kernel Type: Fourth-Order Beta associated (bounded/range)"
  )
})

conditional_kernel_block <- function(object, method = c("summary", "print")) {
  method <- match.arg(method)
  output <- capture.output(if (identical(method, "summary")) {
    summary(object)
  } else {
    print(object)
  })
  trimws(grep(
    "^(Continuous|Unordered Categorical|Ordered Categorical) Kernel Type|^No\\. (Continuous|Unordered Categorical|Ordered Categorical)",
    output, value = TRUE
  ))
}

make_mixed_conditional_bw <- function(kind = c("density", "distribution"),
                                      same = FALSE) {
  kind <- match.arg(kind)
  x <- data.frame(
    xc = seq(0.10, 0.90, length.out = 24L),
    xu = factor(rep(c("a", "b"), 12L)),
    xo = ordered(rep(c("low", "middle", "high"), 8L))
  )
  y <- data.frame(
    yc = seq(0.12, 0.88, length.out = 24L),
    yu = factor(rep(c("u", "v", "w"), 8L)),
    yo = ordered(rep(c("down", "up"), 12L))
  )
  arguments <- list(
    xdat = x, ydat = y, bws = rep(0.2, 6L),
    bandwidth.compute = FALSE, regtype = "lp", degree = 1L,
    uxkertype = if (same && identical(kind, "distribution"))
      "aitchisonaitken" else "liracine",
    oxkertype = if (same) "liracine" else "wangvanryzin",
    oykertype = "liracine"
  )
  if (identical(kind, "density"))
    arguments$uykertype <- if (same) "liracine" else "aitchisonaitken"

  do.call(if (identical(kind, "density")) npcdensbw else npcdistbw,
          arguments)
}

test_that("all conditional kernel families retain explicit X and Y roles", {
  expected <- c(
    "Continuous Kernel Type (Exp. Var.): Second-Order Gaussian",
    "Continuous Kernel Type (Dep. Var.): Second-Order Gaussian",
    "No. Continuous Explanatory Vars.: 1",
    "No. Continuous Dependent Vars.: 1",
    "Unordered Categorical Kernel Type (Exp. Var.): Li and Racine (normalized)",
    "Unordered Categorical Kernel Type (Dep. Var.): Li and Racine (normalized)",
    "No. Unordered Categorical Explanatory Vars.: 1",
    "No. Unordered Categorical Dependent Vars.: 1",
    "Ordered Categorical Kernel Type (Exp. Var.): Li and Racine",
    "Ordered Categorical Kernel Type (Dep. Var.): Li and Racine (normalized)",
    "No. Ordered Categorical Explanatory Vars.: 1",
    "No. Ordered Categorical Dependent Vars.: 1"
  )

  for (kind in c("density", "distribution")) {
    object <- make_mixed_conditional_bw(kind, same = TRUE)
    kind.expected <- expected
    if (identical(kind, "distribution")) {
      kind.expected[5:6] <- paste0(
        c(
          "Unordered Categorical Kernel Type (Exp. Var.): ",
          "Unordered Categorical Kernel Type (Dep. Var.): "
        ),
        "Aitchison and Aitken"
      )
    }
    expect_identical(conditional_kernel_block(object), kind.expected)
    expect_identical(conditional_kernel_block(object, "print"), kind.expected)
  }
})

test_that("conditional kernel sections use one canonical print and summary layout", {
  for (kind in c("density", "distribution")) {
    object <- make_mixed_conditional_bw(kind)
    expect_identical(
      conditional_kernel_block(object, "print"),
      conditional_kernel_block(object, "summary")
    )
  }
})

test_that("the obsolete conditional formatter wrapper is absent", {
  namespace <- environment(genBwKerStrs)
  expect_false(exists("genBwKerStrsXY", envir = namespace, inherits = FALSE))
})

test_that("named categorical matching preserves match.arg semantics", {
  unordered <- c("aitchisonaitken", "liracine")
  ordered <- c("liracine", "wangvanryzin", "racineliyan")

  expect_identical(
    .np_match_arg_named(unordered, unordered, "uxkertype"),
    unordered[[1L]]
  )
  expect_identical(
    .np_match_arg_named("a", unordered, "uykertype"),
    "aitchisonaitken"
  )
  expect_identical(
    .np_match_arg_named("w", ordered, "oxkertype"),
    "wangvanryzin"
  )

  bad.values <- list("invalid", character(), c("liracine", "wangvanryzin"), 1)
  for (value in bad.values) {
    condition <- capture_np_error(
      .np_match_arg_named(value, ordered, "oykertype")
    )
    expect_s3_class(condition, "simpleError")
    expect_match(conditionMessage(condition), "'oykertype'", fixed = TRUE)
    expect_identical(deparse(conditionCall(condition)),
                     "match.arg(oykertype)")
  }
})

test_that("conditional categorical diagnostics name each public argument", {
  x <- data.frame(x = seq(0.05, 0.95, length.out = 24L))
  y <- data.frame(y = seq(0.08, 0.92, length.out = 24L))
  common <- list(
    xdat = x, ydat = y, bws = c(0.2, 0.2),
    bandwidth.compute = FALSE
  )

  for (argument in c("uxkertype", "uykertype", "oxkertype", "oykertype")) {
    call.args <- common
    call.args[[argument]] <- "invalid"
    condition <- capture_np_error(do.call(npcdensbw, call.args))
    expect_s3_class(condition, "simpleError")
    expect_match(conditionMessage(condition),
                 paste0("'", argument, "'"), fixed = TRUE)
    expect_identical(
      deparse(conditionCall(condition)),
      paste0("match.arg(", argument, ")")
    )
  }

  for (argument in c("uxkertype", "oxkertype", "oykertype")) {
    call.args <- common
    call.args[[argument]] <- "invalid"
    condition <- capture_np_error(do.call(npcdistbw, call.args))
    expect_s3_class(condition, "simpleError")
    expect_match(conditionMessage(condition),
                 paste0("'", argument, "'"), fixed = TRUE)
  }
})
