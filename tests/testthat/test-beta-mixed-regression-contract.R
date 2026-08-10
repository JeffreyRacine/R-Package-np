beta_mixed_regression_kernel_args <- function(bws, training, evaluation) {
  bandwidth <- bws[["bw", exact = TRUE]]
  categorical <- bws[["iuno", exact = TRUE]] |
    bws[["iord", exact = TRUE]]

  if (isTRUE(bws[["scaling", exact = TRUE]])) {
    realized <- bws[["bandwidth", exact = TRUE]]
    stopifnot(is.list(realized), length(realized) == 1L)
    bandwidth[categorical] <- as.numeric(realized[[1L]])[categorical]
  }
  arguments <- list(
    bws = bandwidth,
    txdat = training,
    exdat = evaluation,
    bwtype = bws[["type", exact = TRUE]],
    ckertype = "beta",
    ckerorder = bws[["ckerorder", exact = TRUE]],
    ckerbound = bws[["ckerbound", exact = TRUE]],
    ukertype = bws[["ukertype", exact = TRUE]],
    okertype = bws[["okertype", exact = TRUE]],
    return.kernel.weights = TRUE
  )
  if (identical(bws[["ckerbound", exact = TRUE]], "fixed")) {
    arguments$ckerlb <- bws[["ckerlb", exact = TRUE]][
      bws[["icon", exact = TRUE]]
    ]
    arguments$ckerub <- bws[["ckerub", exact = TRUE]][
      bws[["icon", exact = TRUE]]
    ]
  }
  arguments
}

beta_mixed_regression_moments <- function(weights, response) {
  denominator <- colSums(weights)
  normalized <- sweep(weights, 2L, denominator, "/")
  mean <- colSums(normalized * response)
  centered <- response - matrix(
    mean, nrow = nrow(weights), ncol = ncol(weights), byrow = TRUE
  )
  variance <- pmax(colSums(normalized * centered^2), 0)

  list(
    mean = mean,
    stderr = sqrt(variance * colSums(normalized^2)),
    denominator = denominator
  )
}

locate_beta_mixed_regression_sources <- function() {
  roots <- unique(c(
    testthat::test_path("..", ".."),
    testthat::test_path("..", "..", ".."),
    Sys.getenv("R_PACKAGE_DIR", ""),
    Sys.getenv("R_PACKAGE_SOURCE", ""),
    getwd(), file.path(getwd(), "..")
  ))
  roots <- roots[nzchar(roots)]
  roots <- roots[file.exists(file.path(roots, "src", "jksum.c"))]
  if (!length(roots))
    return(NULL)
  roots[[1L]]
}

test_that("mixed beta scalar regression matches canonical kernel-weight oracles", {
  old <- options(np.messages = FALSE, np.categorical.compress = FALSE)
  on.exit(options(old), add = TRUE)

  x <- seq(0.025, 0.975, length.out = 24L)
  training <- data.frame(
    x = x,
    u = factor(rep(letters[1:3], length.out = length(x))),
    o = ordered(rep(1:4, each = 6L))
  )
  response <- sin(4 * x) + 0.12 * as.integer(training$u) +
    0.07 * as.integer(training$o)
  evaluation <- training[c(2L, 5L, 9L, 14L, 20L, 23L), , drop = FALSE]

  for (order in c(2L, 4L, 6L, 8L)) {
    for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
      bandwidth <- if (identical(bwtype, "fixed")) {
        c(0.15, 0.2, 0.25)
      } else {
        c(7, 0.2, 0.25)
      }
      bws <- npregbw(
        xdat = training, ydat = response,
        bws = bandwidth, bandwidth.compute = FALSE,
        regtype = "lc", bwtype = bwtype,
        ckertype = "beta", ckerorder = order,
        ckerbound = "fixed", ckerlb = 0, ckerub = 1
      )
      kernel_arguments <- beta_mixed_regression_kernel_args(
        bws, training, evaluation
      )
      weights <- do.call(npksum, kernel_arguments)$kw
      expected <- beta_mixed_regression_moments(weights, response)

      options(np.categorical.compress = FALSE)
      dense <- npreg(
        bws = bws, txdat = training, tydat = response,
        exdat = evaluation, gradients = TRUE, se = TRUE
      )
      options(np.categorical.compress = TRUE)
      compressed <- npreg(
        bws = bws, txdat = training, tydat = response,
        exdat = evaluation, gradients = TRUE, se = TRUE
      )

      expect_equal(fitted(dense), expected$mean, tolerance = 3e-10)
      expect_equal(se(dense), expected$stderr, tolerance = 3e-10)
      expect_identical(fitted(compressed), fitted(dense))
      expect_identical(se(compressed), se(dense))
      expect_identical(gradients(compressed), gradients(dense))
      expect_identical(compressed$gerr, dense$gerr)

      derivative_arguments <- kernel_arguments
      derivative_arguments$operator <- c(
        "derivative", "normal", "normal"
      )
      derivative_weights <- do.call(npksum, derivative_arguments)$kw
      expected_continuous_gradient <-
        (colSums(derivative_weights * response) -
           expected$mean * colSums(derivative_weights)) /
        expected$denominator
      expect_equal(
        gradients(dense)[, 1L], expected_continuous_gradient,
        tolerance = 4e-10
      )

      unordered_evaluation <- evaluation
      unordered_evaluation$u <- factor(
        rep(levels(training$u)[1L], nrow(evaluation)),
        levels = levels(training$u)
      )
      unordered_arguments <- kernel_arguments
      unordered_arguments$exdat <- unordered_evaluation
      unordered_expected <- beta_mixed_regression_moments(
        do.call(npksum, unordered_arguments)$kw, response
      )
      expect_equal(
        gradients(dense)[, 2L],
        expected$mean - unordered_expected$mean,
        tolerance = 4e-10
      )
      expect_equal(
        dense$gerr[, 2L],
        sqrt(expected$stderr^2 + unordered_expected$stderr^2),
        tolerance = 4e-10
      )

      ordered_evaluation <- evaluation
      ordered_index <- match(
        as.character(ordered_evaluation$o), levels(training$o)
      )
      alternate_index <- ifelse(
        ordered_index == 1L, 2L, ordered_index - 1L
      )
      ordered_evaluation$o <- ordered(
        levels(training$o)[alternate_index], levels = levels(training$o)
      )
      ordered_arguments <- kernel_arguments
      ordered_arguments$exdat <- ordered_evaluation
      ordered_expected <- beta_mixed_regression_moments(
        do.call(npksum, ordered_arguments)$kw, response
      )
      expect_equal(
        gradients(dense)[, 3L],
        (expected$mean - ordered_expected$mean) *
          ifelse(ordered_index == 1L, -1, 1),
        tolerance = 4e-10
      )
      expect_equal(
        dense$gerr[, 3L],
        sqrt(expected$stderr^2 + ordered_expected$stderr^2),
        tolerance = 4e-10
      )
    }
  }
})

test_that("mixed beta regression preserves public routes and objective search", {
  old <- options(np.messages = FALSE, np.categorical.compress = TRUE)
  on.exit(options(old), add = TRUE)

  training <- data.frame(
    x = seq(0.04, 0.96, length.out = 18L),
    u = factor(rep(letters[1:3], 6L)),
    o = ordered(rep(1:3, each = 6L))
  )
  training$y <- cos(3 * training$x) +
    0.1 * as.integer(training$u) + 0.05 * as.integer(training$o)
  evaluation <- training[c(2L, 7L, 12L, 17L), c("x", "u", "o")]
  common <- list(
    bws = c(0.16, 0.2, 0.25), regtype = "lc",
    ckertype = "beta", ckerorder = 6,
    ckerbound = "fixed", ckerlb = 0, ckerub = 1
  )
  direct <- do.call(npreg, c(list(
    txdat = training[c("x", "u", "o")], tydat = training$y,
    exdat = evaluation, gradients = TRUE, se = TRUE
  ), common))
  formula <- do.call(npreg, c(list(
    y ~ x + u + o, data = training, newdata = evaluation,
    gradients = TRUE, se = TRUE
  ), common))
  bws <- do.call(npregbw, c(list(
    xdat = training[c("x", "u", "o")], ydat = training$y,
    bandwidth.compute = FALSE
  ), common))
  object_fit <- npreg(
    bws = bws, txdat = training[c("x", "u", "o")],
    tydat = training$y, exdat = evaluation, gradients = TRUE, se = TRUE
  )
  direct_helper <- np:::.np_regression_direct(
    bws = bws, txdat = training[c("x", "u", "o")],
    tydat = training$y, exdat = evaluation, gradients = TRUE, se = TRUE
  )
  prediction <- predict(formula, newdata = evaluation, se.fit = TRUE)

  expect_equal(fitted(formula), fitted(direct), tolerance = 3e-10)
  expect_equal(fitted(object_fit), fitted(direct), tolerance = 3e-10)
  expect_equal(gradients(formula), gradients(direct), tolerance = 3e-10)
  expect_equal(gradients(object_fit), gradients(direct), tolerance = 3e-10)
  expect_equal(direct_helper$mean, fitted(direct), tolerance = 3e-10)
  expect_equal(direct_helper$grad, gradients(direct), tolerance = 3e-10)
  expect_equal(as.numeric(prediction$fit), fitted(direct), tolerance = 3e-10)
  expect_equal(as.numeric(prediction$se.fit), se(direct), tolerance = 3e-10)

  selected <- npregbw(
    xdat = training[c("x", "u", "o")], ydat = training$y,
    bws = c(0.16, 0.2, 0.25), bandwidth.compute = FALSE,
    regtype = "lc", ckertype = "beta",
    ckerbound = "fixed", ckerlb = 0, ckerub = 1
  )
  expect_true(is.finite(np:::.npregbw_eval_only(
    training[c("x", "u", "o")], training$y, selected,
    invalid.penalty = "dbmax"
  )$objective[[1L]]))
  options(np.categorical.compress = NA)
  expect_error(
    npreg(bws = bws, txdat = training[c("x", "u", "o")],
          tydat = training$y),
    "np.categorical.compress", fixed = TRUE
  )
  expect_error(
    npregbw(
      xdat = training[c("x", "u", "o")], ydat = training$y,
      regtype = "lc", ckertype = "beta",
      ckerbound = "fixed", ckerlb = 0, ckerub = 1
    ),
    "np.categorical.compress", fixed = TRUE
  )
  expect_silent(npreg(
    bws = c(0.16, 0.2, 0.25),
    txdat = training[c("x", "u", "o")], tydat = training$y,
    regtype = "lc", ckertype = "gaussian"
  ))
})

test_that("mixed beta regression is structurally canonical and linear-memory", {
  root <- locate_beta_mixed_regression_sources()
  skip_if(is.null(root), "package sources unavailable")
  ingress <- paste(
    readLines(file.path(root, "src", "np.c"), warn = FALSE),
    collapse = "\n"
  )
  engine <- paste(
    readLines(file.path(root, "src", "jksum.c"), warn = FALSE),
    collapse = "\n"
  )

  expect_false(grepl(
    "beta regression currently supports continuous predictors only",
    ingress, fixed = TRUE
  ))
  expect_match(ingress, "REG_CATCOMPI", fixed = TRUE)
  expect_match(
    engine, "np_beta_regression_categorical_gradients_validated(",
    fixed = TRUE
  )
  expect_match(
    engine, "np_beta_categorical_factor_context_prepare(", fixed = TRUE
  )
  helper_start <- regexpr(
    "np_beta_regression_categorical_gradients_validated(",
    engine, fixed = TRUE
  )[[1L]]
  helper_end <- regexpr(
    "np_beta_regression_gradient_rows_validated(",
    engine, fixed = TRUE
  )[[1L]]
  expect_gt(helper_start, 0L)
  expect_gt(helper_end, helper_start)
  helper <- substr(engine, helper_start, helper_end - 1L)
  expect_false(grepl(
    "plan->num_train * plan->num_eval", helper, fixed = TRUE
  ))
})

test_that("mixed beta native MADS returns the point it evaluated", {
  skip_if_not_installed("crs")
  old <- options(np.messages = FALSE, np.categorical.compress = TRUE)
  on.exit(options(old), add = TRUE)

  set.seed(2323)
  n <- 81L
  training <- data.frame(
    x1 = rbeta(n, 1.2, 1.7),
    x2 = rbeta(n, 1.5, 1.3),
    u = factor(sample(letters[1:3], n, replace = TRUE)),
    o = ordered(sample(1:4, n, replace = TRUE), levels = 1:4)
  )
  response <- sin(3 * training$x1) + 0.35 * training$x2 +
    0.1 * (training$u == "b") - 0.08 * (training$o == "4") +
    rnorm(n, sd = 0.03)
  bandwidth <- npregbw(
    xdat = training,
    ydat = response,
    regtype = "lp",
    degree = c(1L, 1L),
    basis = "glp",
    bernstein.basis = TRUE,
    bwmethod = "cv.ls",
    bwsolver = "mads",
    ckertype = "beta",
    ckerorder = 4L,
    ckerbound = "fixed",
    ckerlb = c(0, 0),
    ckerub = c(1, 1),
    nmulti = 1L,
    nomad.opts = list(MAX_BB_EVAL = 20L)
  )
  replay <- np:::.npregbw_eval_only(
    training, response, bandwidth, invalid.penalty = "baseline"
  )$objective[[1L]]

  expect_identical(as.numeric(bandwidth$fval[[1L]]), as.numeric(replay))
})
