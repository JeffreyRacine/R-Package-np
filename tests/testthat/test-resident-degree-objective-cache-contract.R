test_that("prepared fixed-degree npreg MPI cache keys repeated bandwidths", {
  skip_if_not_installed("crs")
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE,
                      npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(42)
  n <- 120L
  x <- runif(n)
  y <- sin(2 * pi * x) + 0.3 * x + rnorm(n, sd = 0.18)
  xdat <- data.frame(x = x)
  bw <- np::npregbw(
    xdat = xdat,
    ydat = y,
    regtype = "lp",
    degree = 2L,
    bwmethod = "cv.ls",
    bwtype = "fixed",
    bandwidth.compute = FALSE,
    bws = 0.24
  )
  points <- data.frame(
    bw = c(0.24, 0.24, 0.41, 0.24),
    degree = rep(2L, 4L)
  )

  evaluate <- function(cache) {
    options(np.objective.cache = cache)
    prepared <- npRmpi:::.npregbw_prepared_begin(
      xdat = xdat,
      ydat = y,
      bws = bw,
      start.bw = 0.24,
      comm = 1L,
      broadcast = TRUE
    )
    on.exit(npRmpi:::.npregbw_prepared_end(prepared, broadcast = TRUE),
            add = TRUE)
    vapply(seq_len(nrow(points)), function(i) {
      npRmpi:::.npregbw_prepared_eval(
        prepared = prepared,
        bw = points$bw[i],
        degree = points$degree[i],
        broadcast = TRUE
      )[1L]
    }, numeric(1L))
  }

  cached <- evaluate(TRUE)
  uncached <- evaluate(FALSE)
  expect_equal(cached, uncached, tolerance = 1e-12)
  expect_identical(cached[1L], cached[2L])
  expect_false(isTRUE(all.equal(cached[1L], cached[3L], tolerance = 1e-8)))
  expect_identical(cached[1L], cached[4L])
})

test_that("prepared regression and resident density bound degree-key writes", {
  src <- readLines(npRmpi_test_source_path("src", "np.c"), warn = FALSE)
  text <- paste(src, collapse = "\n")
  expect_match(text, "bwm_nn_cache_configure_for_degree_search\\(BANDWIDTH_den_extern")
  expect_match(
    text,
    "prepared_context->degree_key_len = degree_search ?",
    fixed = TRUE
  )
  expect_match(text, "np_conditional_density_nomad_shadow\\.degree_key_len = degree_key_len")
  expect_match(
    text,
    "i < context->degree_key_len"
  )
  expect_match(
    text,
    "i < np_conditional_density_nomad_shadow\\.degree_key_len"
  )
  expect_match(text, "context->num_var \\+ i \\+ 1")
  expect_match(text, "np_conditional_density_nomad_shadow.num_all_var \\+ i \\+ 1")
  expect_match(
    text,
    paste0(
      "if \\(np_conditional_density_nomad_shadow\\.degree_search \\|\\|",
      "\\s+degree_refresh_needed\\)\\s+",
      "degree_refresh_ok = np_mpi_comm1_all_ok\\(degree_refresh_ok\\)"
    ),
    perl = TRUE
  )
  expect_false(grepl("np_regression_nomad_shadow", text, fixed = TRUE))
})
