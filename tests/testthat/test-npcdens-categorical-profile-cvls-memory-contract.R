test_that("npcdens categorical profile CVLS owns a bounded tri-state workspace", {
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
      "np_conditional_density_cvls_categorical_profile_stream\\("
    ),
    source
  )[[1L]]
  start <- tail(starts[starts > 0L], 1L)
  expect_gt(start, 0L)
  tail <- substring(source, start)
  finish <- regexpr(
    paste0(
      "\\nstatic NPConditionalProfileCvStatus\\n",
      "np_conditional_distribution_cvls_categorical_profile_stream\\("
    ),
    tail
  )[[1L]]
  expect_gt(finish, 0L)
  implementation <- substring(tail, 1L, finish - 1L)

  expect_match(
    source,
    "NP_CONDITIONAL_PROFILE_CV_KERNEL_MAX_BYTES",
    fixed = TRUE
  )
  expect_match(
    implementation,
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
  expect_match(implementation, "x_cache_rows = 1", fixed = TRUE)
  expect_match(implementation, "yn_cache_rows = 1", fixed = TRUE)
  expect_match(implementation, "yc_cache_rows = 1", fixed = TRUE)
  expect_match(implementation, "x_cache_rows - 1", fixed = TRUE)
  expect_match(implementation, "yn_cache_rows - 1", fixed = TRUE)
  expect_match(implementation, "yc_cache_rows - 1", fixed = TRUE)
  expect_match(
    implementation,
    "np_hot_loop_check_interrupt",
    fixed = TRUE
  )
  expect_match(
    implementation,
    "Give it first claim on resident",
    fixed = TRUE
  )
  expect_match(
    source,
    "# define NP_CDENS_PROFILE_GROUP_WIDTH 32",
    fixed = TRUE
  )
  expect_match(
    implementation,
    "x_cache_rows - x_scratch_rows + q",
    fixed = TRUE
  )
  expect_match(
    implementation,
    "yn_cache_rows - yn_scratch_rows + q",
    fixed = TRUE
  )
  expect_match(
    implementation,
    "profile_group_width /= 2",
    fixed = TRUE
  )
  expect_match(
    implementation,
    "if(profile_group_width <= 1)",
    fixed = TRUE
  )
  expect_false(grepl(
    "alloc_vecd(nprof_x*nprof_x)",
    implementation,
    fixed = TRUE
  ))
  expect_false(grepl(
    "alloc_vecd(nprof_y*nprof_y)",
    implementation,
    fixed = TRUE
  ))
  expect_false(grepl("counts_y", implementation, fixed = TRUE))
})

test_that("npcdens categorical profile CVLS remains numerically equivalent", {
  set.seed(2026073014L)
  n <- 768L
  x <- data.frame(
    xu = factor(sample(letters[1:4], n, TRUE)),
    xo = ordered(sample(1:6, n, TRUE))
  )
  y <- data.frame(
    yu = factor(sample(LETTERS[1:3], n, TRUE)),
    yo = ordered(sample(1:5, n, TRUE))
  )
  state <- npcdensbw(
    xdat = x,
    ydat = y,
    bwmethod = "cv.ls",
    bwtype = "fixed",
    regtype = "lc",
    uxkertype = "liracine",
    uykertype = "liracine",
    oxkertype = "racineliyan",
    oykertype = "wangvanryzin",
    bws = c(0.19, 0.24, 0.29, 0.34),
    bandwidth.compute = FALSE,
    itmax = 1L
  )
  old <- options(np.categorical.compress = TRUE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  compressed <- np:::.npcdensbw_eval_only(x, y, state)$objective
  options(np.categorical.compress = FALSE)
  canonical <- np:::.npcdensbw_eval_only(x, y, state)$objective

  expect_equal(compressed, canonical, tolerance = 1e-12)
})
