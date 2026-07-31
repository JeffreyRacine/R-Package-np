test_that("npcdist categorical profile CVLS owns a bounded MPI tri-state workspace", {
  candidates <- c(
    test_path("..", ".."),
    test_path("..", "..", ".."),
    Sys.getenv("R_PACKAGE_SOURCE", unset = "")
  )
  source_root <- candidates[
    file.exists(file.path(candidates, "src", "jksum.c"))
  ][1L]
  skip_if(is.na(source_root), "package source unavailable")

  source <- paste(
    readLines(file.path(source_root, "src", "jksum.c"), warn = FALSE),
    collapse = "\n"
  )
  starts <- gregexpr(
    paste0(
      "static NPConditionalProfileCvStatus\\n",
      "np_conditional_distribution_cvls_categorical_profile_stream\\("
    ),
    source
  )[[1L]]
  start <- tail(starts[starts > 0L], 1L)
  expect_gt(start, 0L)
  tail <- substring(source, start)
  finish <- regexpr(
    "\\n/\\*\\n \\* Return 2 only when no optional rank-owned slab",
    tail
  )[[1L]]
  expect_gt(finish, 0L)
  implementation <- substring(tail, 1L, finish - 1L)

  expect_match(
    source,
    "np_conditional_profile_cv_serial_kernel_max_bytes",
    fixed = TRUE
  )
  expect_match(
    source,
    "np_conditional_profile_cv_mpi_kernel_max_bytes",
    fixed = TRUE
  )
  expect_match(
    source,
    "NP_CONDITIONAL_PROFILE_CV_NOT_APPLICABLE",
    fixed = TRUE
  )
  expect_match(
    source,
    "if(profile_status == NP_CONDITIONAL_PROFILE_CV_FAILURE)",
    fixed = TRUE
  )
  expect_match(implementation, "np_size_mul_checked", fixed = TRUE)
  expect_match(implementation, "np_size_add_checked", fixed = TRUE)
  expect_match(implementation, "np_size_array_bytes_checked", fixed = TRUE)
  expect_match(implementation, "np_native_malloc_array", fixed = TRUE)
  expect_match(
    implementation,
    "np_categorical_profile_tile_fill_prevalidated",
    fixed = TRUE
  )
  expect_match(implementation, "use_x_cache", fixed = TRUE)
  expect_match(implementation, "use_y_cache", fixed = TRUE)
  expect_match(implementation, "x_cache_rows = 1", fixed = TRUE)
  expect_match(implementation, "y_cache_rows = 1", fixed = TRUE)
  expect_match(implementation, "y_cache_rows - 1", fixed = TRUE)
  expect_match(
    implementation,
    "np_hot_loop_check_interrupt",
    fixed = TRUE
  )
  expect_match(
    source,
    "np_conditional_profile_cache_priority_compare",
    fixed = TRUE
  )
  expect_match(
    implementation,
    "MPI_Allreduce(&local_preflight_status",
    fixed = TRUE
  )
  expect_match(
    implementation,
    "NP_RMPI_INJECT_CDIST_PROFILE_FAIL_RANK",
    fixed = TRUE
  )
  expect_false(grepl(
    "alloc_vecd(nprof_x*nprof_x)",
    implementation,
    fixed = TRUE
  ))
  expect_false(grepl(
    "alloc_vecd(nprof_ey*nprof_ty)",
    implementation,
    fixed = TRUE
  ))
  expect_false(grepl("counts_ty", implementation, fixed = TRUE))
})

test_that("npcdist categorical profile CVLS remains numerically equivalent", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(), add = TRUE)
  set.seed(2026073007L)
  n <- 768L
  x <- data.frame(
    xu = factor(sample(letters[1:4], n, TRUE)),
    xo = ordered(sample(1:6, n, TRUE))
  )
  y <- data.frame(
    yu = factor(sample(LETTERS[1:3], n, TRUE)),
    yo = ordered(sample(1:5, n, TRUE))
  )
  state <- npcdistbw(
    xdat = x,
    ydat = y,
    bwmethod = "cv.ls",
    bwtype = "fixed",
    regtype = "lc",
    uxkertype = "liracine",
    oxkertype = "racineliyan",
    oykertype = "wangvanryzin",
    bws = c(0.19, 0.24, 0.29, 0.34),
    bandwidth.compute = FALSE,
    itmax = 1L
  )
  old <- options(np.categorical.compress = TRUE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  compressed <- npRmpi:::.npcdistbw_eval_only(
    xdat = x,
    ydat = y,
    bws = state,
    do.full.integral = TRUE
  )$objective
  options(np.categorical.compress = FALSE)
  canonical <- npRmpi:::.npcdistbw_eval_only(
    xdat = x,
    ydat = y,
    bws = state,
    do.full.integral = TRUE
  )$objective

  expect_equal(compressed, canonical, tolerance = 1e-12)
})
