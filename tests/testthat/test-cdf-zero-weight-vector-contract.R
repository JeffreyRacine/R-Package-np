locate_cdf_vector_source <- function() {
  candidates <- c(
    test_path("..", "..", "src", "jksum.c"),
    test_path("..", "..", "..", "src", "jksum.c"),
    file.path(Sys.getenv("R_PACKAGE_DIR", ""), "src", "jksum.c"),
    file.path(Sys.getenv("R_PACKAGE_SOURCE", ""), "src", "jksum.c"),
    file.path(getwd(), "src", "jksum.c"),
    file.path(getwd(), "..", "src", "jksum.c")
  )
  candidates <- unique(candidates[nzchar(candidates)])
  hits <- candidates[file.exists(candidates)]
  if (!length(hits)) NULL else hits[[1L]]
}

test_that("ordinary Gaussian CDF weighted rows preserve zero-lane semantics", {
  src_file <- locate_cdf_vector_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  source <- paste(readLines(src_file, warn = FALSE), collapse = "\n")
  start <- regexpr(
    "static int NP_NOINLINE np_ckernelv_accel_try(",
    source,
    fixed = TRUE
  )
  stop <- regexpr("void np_ckernelv(", source, fixed = TRUE)
  expect_gt(start, 0L)
  expect_gt(stop, start)
  helper <- substr(source, start, stop - 1L)

  expect_match(helper, "(KERNEL != 30 || num_xt <= 8192)", fixed = TRUE)
  expect_match(helper, "if(KERNEL == 30)", fixed = TRUE)
  expect_match(helper, "xw[i] != 0.0 && xw[i] != 1.0", fixed = TRUE)
  expect_match(helper, "if(!have_nonunit_weight)", fixed = TRUE)
  expect_match(helper, "np_accel_gauss_val[i] = 1.0", fixed = TRUE)
  expect_match(
    helper,
    "np_accel_gauss_val, 1, xw, 1, result, 1",
    fixed = TRUE
  )
  expect_match(
    helper,
    "if(!np_accel_gauss_has_zero_weight(xw, num_xt))",
    fixed = TRUE
  )
  expect_lt(
    regexpr("if(KERNEL == 30)", helper, fixed = TRUE),
    regexpr(
      "if(!np_accel_gauss_has_zero_weight(xw, num_xt))",
      helper,
      fixed = TRUE
    )
  )
  expect_false(grepl("malloc(", helper, fixed = TRUE))
  expect_false(grepl("calloc(", helper, fixed = TRUE))
  expect_false(grepl("realloc(", helper, fixed = TRUE))
})

test_that("accelerated Gaussian CVCDF remains numerically tied to scalar CVCDF", {
  owner <- environmentName(environment(npudistbw))
  skip_if(
    identical(owner, "npRmpi"),
    "active-pool numerical parity is covered by the MPI qualification sentinel"
  )

  set.seed(20260727L)
  dat <- data.frame(x1 = rnorm(300L), x2 = rnorm(300L))
  bw <- npudistbw(
    dat = dat,
    bws = c(0.35, 0.38),
    bandwidth.compute = FALSE,
    bwmethod = "cv.cdf"
  )
  evaluate <- function(accelerate) {
    old <- options(np.macMseries.accelerate = accelerate)
    on.exit(options(old), add = TRUE)
    getFromNamespace("npudistbw.dbandwidth", owner)(
      dat = dat,
      bws = bw,
      bandwidth.compute = TRUE,
      eval.only = TRUE,
      nmulti = 1L,
      powell.remin = FALSE
    )$fval[[1L]]
  }

  scalar <- evaluate(FALSE)
  accelerated <- evaluate(TRUE)
  expect_true(is.finite(scalar))
  expect_true(is.finite(accelerated))
  expect_equal(accelerated, scalar, tolerance = 1e-10)
})
