high_occupancy_fixture <- function() {
  set.seed(20260815L)
  n <- 40L
  x <- data.frame(
    x1 = sort(runif(n, -1, 1)),
    x2 = runif(n, -0.2, 0.2)
  )
  y <- 0.7 + 0.5 * x$x1 - 0.3 * x$x2 + 0.25 * x$x1 * x$x2 +
    0.08 * x$x1^2 + rnorm(n, sd = 0.04)
  distances <- sort(outer(
    x$x1, x$x1, function(a, b) abs(a - b)
  )[lower.tri(diag(n))])
  total <- n * (n - 1L) / 2L
  threshold <- ceiling(9 * total / 10)
  radius <- function(count) mean(distances[c(count, count + 1L)])
  list(
    x = x,
    y = y,
    n = n,
    total = total,
    threshold = threshold,
    radius = radius,
    full2 = diff(range(x$x2)) * 1.01
  )
}

high_occupancy_eval <- function(fixture, degree, bandwidth,
                                method = "cv.ls", objective = "ls",
                                kernel = "uniform", order = 2L,
                                bernstein = TRUE,
                                x = fixture$x) {
  arguments <- list(
    xdat = x,
    ydat = fixture$y,
    bws = bandwidth,
    bandwidth.compute = FALSE,
    bwtype = "fixed",
    bwmethod = method,
    ckertype = kernel,
    regtype = "lp",
    degree = degree,
    bernstein.basis = bernstein
  )
  if (identical(kernel, "epanechnikov"))
    arguments$ckerorder <- order
  bw <- do.call(npregbw, arguments)
  core <- getFromNamespace(".npregbw_call_fixed_degree_core", "np")

  options(np.tree = TRUE)
  Sys.setenv(NP_LP_TREE_ORACLE = "1")
  messages <- capture.output(
    tree <- core(xdat = x, ydat = fixture$y, bws = bw,
                 eval.only = TRUE, objective = objective),
    type = "message"
  )
  Sys.unsetenv("NP_LP_TREE_ORACLE")
  options(np.tree = FALSE)
  dense <- core(xdat = x, ydat = fixture$y, bws = bw,
                eval.only = TRUE, objective = objective)
  list(
    sparse = any(grepl("NP_LP_TREE_ORACLE", messages, fixed = TRUE)),
    tree = tree,
    dense = dense,
    messages = messages
  )
}

test_that("wide fixed CVLS transfers only at the exact occupancy boundary", {
  old <- options(np.messages = FALSE, np.macMseries.accelerate = FALSE)
  on.exit({
    options(old)
    Sys.unsetenv("NP_LP_TREE_ORACLE")
  }, add = TRUE)
  fixture <- high_occupancy_fixture()

  cells <- list(
    below = list(
      degree = c(3L, 3L),
      bandwidth = c(fixture$radius(fixture$threshold - 1L), fixture$full2),
      sparse = TRUE
    ),
    exact = list(
      degree = c(3L, 3L),
      bandwidth = c(fixture$radius(fixture$threshold), fixture$full2),
      sparse = FALSE
    ),
    above = list(
      degree = c(3L, 3L),
      bandwidth = c(fixture$radius(fixture$threshold + 1L), fixture$full2),
      sparse = FALSE
    ),
    lower_width = list(
      degree = c(3L, 2L),
      bandwidth = c(fixture$radius(fixture$threshold), fixture$full2),
      sparse = TRUE
    ),
    two_partial = list(
      degree = c(3L, 3L),
      bandwidth = c(fixture$radius(fixture$threshold), fixture$full2 * 0.6),
      sparse = TRUE
    ),
    full = list(
      degree = c(3L, 3L),
      bandwidth = c(diff(range(fixture$x$x1)) * 1.01, fixture$full2),
      sparse = FALSE
    )
  )

  for (label in names(cells)) {
    cell <- cells[[label]]
    result <- high_occupancy_eval(
      fixture, cell$degree, cell$bandwidth, kernel = "uniform"
    )
    expect_identical(result$sparse, cell$sparse, info = label)
    expect_equal(result$tree$objective, result$dense$objective,
                 tolerance = 2e-11, info = label)
    expect_identical(result$tree$num.feval, result$dense$num.feval,
                     info = label)
  }

  epanechnikov <- high_occupancy_eval(
    fixture,
    c(3L, 3L),
    c(fixture$radius(fixture$threshold), fixture$full2) / sqrt(5),
    kernel = "epanechnikov"
  )
  expect_false(epanechnikov$sparse)
  expect_equal(epanechnikov$tree$objective, epanechnikov$dense$objective,
               tolerance = 2e-11)

  raw <- high_occupancy_eval(
    fixture,
    c(3L, 3L),
    c(fixture$radius(fixture$threshold), fixture$full2),
    kernel = "uniform",
    bernstein = FALSE
  )
  expect_false(raw$sparse)
  expect_equal(raw$tree$objective, raw$dense$objective, tolerance = 2e-11)
})

test_that("high-occupancy transfer preserves objective and kernel exclusions", {
  old <- options(np.messages = FALSE, np.macMseries.accelerate = FALSE)
  on.exit({
    options(old)
    Sys.unsetenv("NP_LP_TREE_ORACLE")
  }, add = TRUE)
  fixture <- high_occupancy_fixture()
  bandwidth <- c(fixture$radius(fixture$threshold), fixture$full2)

  cases <- list(
    cvaic = list(method = "cv.aic", objective = "ls",
                 kernel = "uniform", order = 2L),
    cvks = list(method = "cv.ls", objective = "ks",
                kernel = "uniform", order = 2L),
    epanechnikov4 = list(method = "cv.ls", objective = "ls",
                         kernel = "epanechnikov", order = 4L)
  )
  for (label in names(cases)) {
    case <- cases[[label]]
    result <- high_occupancy_eval(
      fixture,
      c(3L, 3L),
      bandwidth / if (case$kernel == "epanechnikov") sqrt(9) else 1,
      method = case$method,
      objective = case$objective,
      kernel = case$kernel,
      order = case$order
    )
    expect_true(result$sparse, info = label)
    expect_equal(result$tree$objective, result$dense$objective,
                 tolerance = 2e-11, info = label)
    expect_identical(result$tree$num.feval, result$dense$num.feval,
                     info = label)
  }

  lsq.arguments <- list(
    xdat = fixture$x,
    ydat = fixture$y,
    scale = rep.int(1, fixture$n),
    tau = 0.5,
    bws = bandwidth,
    regtype = "lp",
    degree = c(3L, 3L),
    bernstein.basis = TRUE,
    ckertype = "uniform",
    bandwidth.compute = FALSE
  )
  options(np.tree = TRUE)
  Sys.setenv(NP_LP_TREE_ORACLE = "1")
  lsq.messages <- capture.output(
    lsq.tree <- do.call(nplsqregbw, lsq.arguments),
    type = "message"
  )
  Sys.unsetenv("NP_LP_TREE_ORACLE")
  options(np.tree = FALSE)
  lsq.dense <- do.call(nplsqregbw, lsq.arguments)
  expect_true(any(grepl("NP_LP_TREE_ORACLE", lsq.messages, fixed = TRUE)))
  expect_equal(lsq.tree$objective, lsq.dense$objective, tolerance = 2e-11)
  expect_identical(lsq.tree$num.feval.fast, lsq.dense$num.feval.fast)

  mixed.x <- transform(
    fixture$x,
    group = factor(rep(c("a", "b"), length.out = fixture$n))
  )
  mixed <- high_occupancy_eval(
    fixture,
    c(3L, 3L),
    c(bandwidth, 0.25),
    kernel = "uniform",
    x = mixed.x
  )
  expect_true(mixed$sparse)
  expect_equal(mixed$tree$objective, mixed$dense$objective,
               tolerance = 2e-11)
})

test_that("high-occupancy certificate borrows and restores bounded workspace", {
  roots <- unique(c(
    testthat::test_path("..", ".."),
    testthat::test_path("..", "..", ".."),
    Sys.getenv("R_PACKAGE_DIR", ""),
    Sys.getenv("R_PACKAGE_SOURCE", ""),
    getwd(),
    file.path(getwd(), "..")
  ))
  paths <- file.path(roots[nzchar(roots)], "src", "jksum.c")
  paths <- paths[file.exists(paths)]
  skip_if(!length(paths), "package sources unavailable")
  lines <- readLines(paths[[1L]], warn = FALSE)
  source <- paste(lines, collapse = "\n")

  expect_match(source, "NP_REG_CV_LP_DENSE_MIN_TERMS = 10", fixed = TRUE)
  expect_match(source, "NP_REG_CV_LP_DENSE_OCCUPANCY_NUMERATOR = 9",
               fixed = TRUE)
  expect_match(source, "NP_REG_CV_LP_DENSE_OCCUPANCY_DENOMINATOR = 10",
               fixed = TRUE)
  expect_match(source, "partial_count > 1", fixed = TRUE)
  expect_match(source, "if(terms[l] != 0)", fixed = TRUE)
  expect_match(source, "sorted = basis[0];", fixed = TRUE)
  expect_match(
    source,
    "sorted[i] = matrix_X_continuous[partial_dimension][i];",
    fixed = TRUE
  )
  expect_match(source, "sorted[i] = 1.0;", fixed = TRUE)
  expect_match(source, "F77_CALL(dlasrt)", fixed = TRUE)
  expect_match(source,
               "certificate_ok = np_reg_cv_support_sort_increasing(sorted, num_obs);",
               fixed = TRUE)
  expect_match(source, "single exit", fixed = TRUE)
  helper <- np_test_extract_c_function(
    lines, "np_reg_fixed_tree_dense_high_occupancy_admitted"
  )
  expect_false(grepl("malloc(", helper, fixed = TRUE))
  expect_false(grepl("calloc(", helper, fixed = TRUE))
  expect_false(grepl("realloc(", helper, fixed = TRUE))
  expect_false(grepl("MPI_", helper, fixed = TRUE))
  expect_false(grepl("R_CheckUserInterrupt", helper, fixed = TRUE))
  expect_false(grepl("num_obs*(size_t)num_obs", helper, fixed = TRUE))
  expect_false(grepl("qsort(sorted", source, fixed = TRUE))
  expect_false(grepl("NP_REG_CV_SUPPORT_ORDER_", source, fixed = TRUE))
  expect_false(grepl("A10", source, fixed = TRUE))
})
