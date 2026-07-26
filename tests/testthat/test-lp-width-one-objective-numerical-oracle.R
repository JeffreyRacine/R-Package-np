.np_width_one_manual_gaussian_objective <- function(xdat,
                                                    ydat,
                                                    bws,
                                                    bwmethod) {
  x <- as.matrix(xdat)
  n <- nrow(x)
  log_weights <- matrix(0.0, nrow = n, ncol = n)

  for (j in seq_len(ncol(x))) {
    scaled_distance <- outer(x[, j], x[, j], "-") / bws[j]
    log_weights <- log_weights - 0.5 * scaled_distance^2
  }

  weights <- exp(log_weights)
  if (identical(bwmethod, "cv.ls"))
    diag(weights) <- 0.0

  denominator <- rowSums(weights)
  fitted <- drop(weights %*% ydat) / denominator
  mean_squared_error <- mean((ydat - fitted)^2)

  if (identical(bwmethod, "cv.ls"))
    return(mean_squared_error)

  trace_hat <- sum(diag(weights) / denominator)
  log(mean_squared_error) +
    (1.0 + trace_hat / n) / (1.0 - (trace_hat + 2.0) / n)
}

test_that("width-one regression objectives match an independent scalar oracle", {
  package <- if ("npRmpi" %in% loadedNamespaces()) "npRmpi" else "np"
  npregbw_owner <- getFromNamespace("npregbw", package)
  evaluate_objective <- getFromNamespace(".npregbw_eval_only", package)
  local_regression <- if (identical(package, "npRmpi"))
    getFromNamespace(".npRmpi_with_local_regression", package) else NULL
  run_local <- function(fun) {
    if (is.null(local_regression))
      fun()
    else
      local_regression(fun())
  }

  old_options <- options(np.tree = FALSE, np.messages = FALSE)
  on.exit(options(old_options), add = TRUE)

  set.seed(2026072506L)
  n <- 67L
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(2 * pi * x1) * cos(2 * pi * x2) +
    0.15 * x1 + rnorm(n, sd = 0.12)
  fixtures <- list(
    list(x = data.frame(x1 = x1), bws = 0.19),
    list(x = data.frame(x1 = x1, x2 = x2), bws = c(0.22, 0.27))
  )

  for (fixture in fixtures) {
    for (bwmethod in c("cv.ls", "cv.aic")) {
      make_bw <- function(regtype, bernstein = FALSE) {
        args <- list(
          xdat = fixture$x,
          ydat = y,
          regtype = regtype,
          bwmethod = bwmethod,
          bwtype = "fixed",
          ckertype = "gaussian",
          ckerorder = 2L,
          bws = fixture$bws,
          bandwidth.compute = FALSE
        )
        if (identical(regtype, "lp")) {
          args$degree <- rep.int(0L, ncol(fixture$x))
          args$degree.select <- "manual"
          args$basis <- "glp"
          args$bernstein.basis <- bernstein
        }
        run_local(function() do.call(npregbw_owner, args))
      }

      bandwidths <- list(
        lc = make_bw("lc"),
        raw_lp0 = make_bw("lp", FALSE),
        bernstein_lp0 = make_bw("lp", TRUE)
      )
      native <- vapply(
        bandwidths,
        function(bw) run_local(
          function() evaluate_objective(fixture$x, y, bw)$objective
        ),
        numeric(1L)
      )
      reference <- .np_width_one_manual_gaussian_objective(
        fixture$x,
        y,
        fixture$bws,
        bwmethod
      )
      scale <- max(1.0, abs(reference))

      expect_identical(unname(native["lc"]), unname(native["raw_lp0"]))
      expect_identical(
        unname(native["raw_lp0"]),
        unname(native["bernstein_lp0"])
      )
      expect_lt(
        max(abs(native - reference)),
        128.0 * .Machine$double.eps * scale,
        label = paste(
          bwmethod,
          ncol(fixture$x),
          "continuous predictor(s)"
        )
      )
    }
  }
})
