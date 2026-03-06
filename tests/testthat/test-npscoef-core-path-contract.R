test_that("npscoef core path does not reference full-matrix helpers", {
  fn.body <- paste(deparse(body(np:::npscoef.scbandwidth)), collapse = "\n")
  expect_no_match(fn.body, "\\.npscoef_weight_matrix\\(")
  expect_no_match(fn.body, "\\.np_kernel_weights_direct")
  expect_no_match(fn.body, "npscoefhat\\(")
})

test_that("ordinary npscoef fits do not hit smoothcoef full-matrix helper", {
  old.opt <- getOption("np.messages")
  options(np.messages = FALSE)
  on.exit(options(np.messages = old.opt), add = TRUE)

  set.seed(20260306)
  n <- 70L
  x <- runif(n)
  z <- runif(n, min = -1.5, max = 1.5)
  y <- (0.3 + sin(1.4 * z)) * x + 0.15 * cos(z) + rnorm(n, sd = 0.04)

  assign(".np_trace_hits_npscoef", 0L, envir = .GlobalEnv)
  on.exit(rm(".np_trace_hits_npscoef", envir = .GlobalEnv), add = TRUE)
  trace(".npscoef_weight_matrix",
        tracer = quote(assign(".np_trace_hits_npscoef",
                              get(".np_trace_hits_npscoef", envir = .GlobalEnv) + 1L,
                              envir = .GlobalEnv)),
        print = FALSE,
        where = asNamespace("np"))
  on.exit(untrace(".npscoef_weight_matrix", where = asNamespace("np")), add = TRUE)

  run_case <- function(regtype) {
    bw.args <- list(
      xdat = x,
      zdat = z,
      ydat = y,
      bws = 0.28,
      regtype = regtype,
      bandwidth.compute = FALSE
    )
    if (identical(regtype, "lp")) {
      bw.args$basis <- "glp"
      bw.args$degree <- 2L
    }
    bw <- do.call(npscoefbw, bw.args)
    assign(".np_trace_hits_npscoef", 0L, envir = .GlobalEnv)
    fit <- npscoef(
      bws = bw,
      txdat = data.frame(x = x),
      tydat = y,
      tzdat = data.frame(z = z),
      gradients = TRUE,
      errors = TRUE,
      iterate = FALSE,
      betas = TRUE
    )
    expect_s3_class(fit, "smoothcoefficient")
    expect_identical(get(".np_trace_hits_npscoef", envir = .GlobalEnv), 0L, info = regtype)
  }

  run_case("lc")
  run_case("ll")
  run_case("lp")
})
