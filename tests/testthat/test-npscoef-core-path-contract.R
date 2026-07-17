test_that("npscoef core path does not reference full-matrix helpers", {
  fn.body <- paste(deparse(body(np:::npscoef.scbandwidth)), collapse = "\n")
  expect_no_match(fn.body, "\\.npscoef_weight_matrix\\(")
  expect_no_match(fn.body, "\\.np_kernel_weights_direct")
  expect_no_match(fn.body, "npscoefhat\\(")
})

test_that("npscoefhat apply matches matrix output under bounded chunks", {
  old.opt <- getOption("np.messages")
  old.chunk <- getOption("np.scoef.weight.chunk.size")
  options(np.messages = FALSE, np.scoef.weight.chunk.size = 7L)
  on.exit({
    options(np.messages = old.opt)
    options(np.scoef.weight.chunk.size = old.chunk)
  }, add = TRUE)

  set.seed(20260306)
  n <- 70L
  x <- runif(n)
  z <- runif(n, min = -1.5, max = 1.5)
  y <- (0.3 + sin(1.4 * z)) * x + 0.15 * cos(z) + rnorm(n, sd = 0.04)

  assign(".np_trace_hits_npscoef_chunk", 0L, envir = .GlobalEnv)
  assign(".np_trace_sizes_npscoef_chunk", integer(), envir = .GlobalEnv)
  on.exit(rm(list = c(".np_trace_hits_npscoef_chunk", ".np_trace_sizes_npscoef_chunk"),
             envir = .GlobalEnv), add = TRUE)
  trace(".npscoef_effective_weight_chunk",
        tracer = quote({
          assign(".np_trace_hits_npscoef_chunk",
                 get(".np_trace_hits_npscoef_chunk", envir = .GlobalEnv) + 1L,
                 envir = .GlobalEnv)
          assign(".np_trace_sizes_npscoef_chunk",
                 c(get(".np_trace_sizes_npscoef_chunk", envir = .GlobalEnv),
                   length(eval.indices)),
                 envir = .GlobalEnv)
        }),
        print = FALSE,
        where = asNamespace("np"))
  on.exit(untrace(".npscoef_effective_weight_chunk", where = asNamespace("np")), add = TRUE)

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
    assign(".np_trace_hits_npscoef_chunk", 0L, envir = .GlobalEnv)
    assign(".np_trace_sizes_npscoef_chunk", integer(), envir = .GlobalEnv)
    common.args <- list(
      bws = bw,
      txdat = data.frame(x = x),
      tzdat = data.frame(z = z),
      exdat = data.frame(x = x),
      ezdat = data.frame(z = z)
    )
    H <- do.call(npscoefhat, c(common.args, list(output = "matrix")))
    applied <- do.call(npscoefhat, c(common.args, list(y = y, output = "apply")))
    expect_equal(as.vector(H %*% y), as.vector(applied), tolerance = 1e-11)

    chunk.sizes <- get(".np_trace_sizes_npscoef_chunk", envir = .GlobalEnv)
    if (identical(regtype, "lc")) {
      expect_gt(get(".np_trace_hits_npscoef_chunk", envir = .GlobalEnv), 0L)
      expect_true(length(chunk.sizes) > 0L)
      expect_true(all(chunk.sizes >= 1L & chunk.sizes <= 7L))
    }
  }

  run_case("lc")
  run_case("ll")
  run_case("lp")
})
