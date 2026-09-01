locate_lp_source <- function(name) {
  candidates <- c(
    test_path("..", "..", "src", name),
    test_path("..", "..", "..", "src", name),
    file.path(Sys.getenv("R_PACKAGE_DIR", ""), "src", name),
    file.path(Sys.getenv("R_PACKAGE_SOURCE", ""), "src", name),
    file.path(getwd(), "src", name),
    file.path(getwd(), "..", "src", name)
  )
  candidates <- unique(candidates[nzchar(candidates)])
  hits <- candidates[file.exists(candidates)]
  if (length(hits) == 0L) NULL else hits[[1L]]
}

test_that("general LP fit solves only requested variance directions", {
  source_file <- locate_lp_source("jksum.c")
  skip_if(is.null(source_file), "LP source unavailable in this test context")
  source <- paste(readLines(source_file, warn = FALSE), collapse = "\n")

  expect_true(grepl("int variance_nrhs = 1;", source, fixed = TRUE))
  expect_true(grepl("np_glp_gradient_direction_active", source, fixed = TRUE))
  expect_true(grepl("np_lp_variance_quadratic", source, fixed = TRUE))
  expect_false(grepl(
    "rhs_source[i+l*owner->nterms] =\n          i == l ? 1.0 : 0.0",
    source,
    fixed = TRUE
  ))
})

test_that("directional LP variances agree with independent fixed-bandwidth WLS", {
  set.seed(271828)
  n <- 83L
  training <- data.frame(
    x1 = runif(n, 0.06, 0.94),
    x2 = runif(n, 0.08, 0.92)
  )
  response <- sin(2 * pi * training$x1) * cos(2 * pi * training$x2) +
    rnorm(n, sd = 0.08)
  evaluation <- data.frame(
    x1 = seq(0.16, 0.84, length.out = 5L),
    x2 = seq(0.78, 0.22, length.out = 5L)
  )
  degree <- c(2L, 1L)
  w_lp <- get("W.lp", envir = environment(npreg), inherits = FALSE)

  for (kernel in c("gaussian", "beta")) {
    for (bernstein in c(FALSE, TRUE)) {
      common <- list(
        bws = c(0.24, 0.27),
        txdat = training,
        exdat = evaluation,
        bwtype = "fixed",
        ckertype = kernel,
        ckerorder = if (identical(kernel, "beta")) 8L else 2L
      )
      if (identical(kernel, "beta")) {
        common$ckerbound <- "fixed"
        common$ckerlb <- c(0, 0)
        common$ckerub <- c(1, 1)
      }
      fit <- do.call(npreg, c(
        common,
        list(
          tydat = response,
          gradients = TRUE,
          se = TRUE,
          regtype = "lp",
          degree = degree,
          basis = "glp",
          bernstein.basis = bernstein
        )
      ))
      weights <- do.call(npksum, c(
        common,
        list(return.kernel.weights = TRUE)
      ))$kw
      if (!is.matrix(weights))
        weights <- matrix(weights, nrow = n, ncol = nrow(evaluation))

      design <- w_lp(
        xdat = training, degree = degree, basis = "glp",
        Bernstein = bernstein
      )
      evaluation_design <- w_lp(
        xdat = training, exdat = evaluation, degree = degree, basis = "glp",
        Bernstein = bernstein
      )
      derivatives <- lapply(seq_along(degree), function(coordinate) {
        order <- integer(length(degree))
        order[[coordinate]] <- 1L
        w_lp(
          xdat = training, exdat = evaluation, degree = degree,
          gradient.vec = order, basis = "glp", Bernstein = bernstein
        )
      })
      training.common <- common
      training.common$exdat <- NULL
      training.fit <- do.call(npreg, c(
        training.common,
        list(
          tydat = response,
          gradients = FALSE,
          se = FALSE,
          regtype = "lp",
          degree = degree,
          basis = "glp",
          bernstein.basis = bernstein
        )
      ))
      residual <- response - fitted(training.fit)

      expected_se <- numeric(nrow(evaluation))
      expected_gerr <- matrix(0, nrow(evaluation), length(degree))
      for (row in seq_len(nrow(evaluation))) {
        weight <- weights[, row]
        gram <- crossprod(design, design * weight)
        projection <- solve(gram, evaluation_design[row, ])
        influence <- weight * drop(design %*% projection)
        expected_se[[row]] <- sqrt(sum((influence * residual)^2))
        for (coordinate in seq_along(degree)) {
          derivative_projection <- solve(
            gram,
            derivatives[[coordinate]][row, ]
          )
          derivative_influence <-
            weight * drop(design %*% derivative_projection)
          expected_gerr[row, coordinate] <-
            sqrt(sum((derivative_influence * residual)^2))
        }
      }

      label <- paste(kernel, if (bernstein) "Bernstein" else "raw")
      expect_equal(unname(se(fit)), expected_se, tolerance = 2e-8,
                   info = label)
      expect_equal(
        unname(gradients(fit, se = TRUE)),
        expected_gerr,
        tolerance = 2e-7,
        info = label
      )
    }
  }
})
