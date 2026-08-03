test_that("continuous scalar kernel logs have one family-neutral owner", {
  candidates <- unique(c(
    normalizePath(file.path("..", ".."), mustWork = FALSE),
    normalizePath(".", mustWork = FALSE)
  ))
  root <- candidates[vapply(candidates, function(path) {
    file.exists(file.path(path, "src", "kernel_registry.c")) &&
      file.exists(file.path(path, "src", "jksum.c"))
  }, logical(1))][1L]
  if (is.na(root)) root <- ""
  skip_if_not(file.exists(file.path(root, "src", "kernel_registry.c")))

  registry.implementation <- paste(
    readLines(file.path(root, "src", "kernel_registry.c"), warn = FALSE),
    collapse = "\n"
  )
  registry.header <- paste(
    readLines(file.path(root, "src", "kernel_registry.h"), warn = FALSE),
    collapse = "\n"
  )
  registry <- paste(registry.implementation, registry.header, sep = "\n")
  engine <- paste(
    readLines(file.path(root, "src", "jksum.c"), warn = FALSE),
    collapse = "\n"
  )
  native.files <- list.files(
    file.path(root, "src"), pattern = "\\.[ch]$", full.names = TRUE
  )
  native <- paste(vapply(native.files, function(path) {
    paste(readLines(path, warn = FALSE), collapse = "\n")
  }, character(1)), collapse = "\n")

  expect_identical(
    lengths(regmatches(
      registry,
      gregexpr("np_continuous_kernel_scalar_log(", registry, fixed = TRUE)
    )),
    2L
  )
  expect_match(
    engine, "np_continuous_kernel_scalar_log(", fixed = TRUE
  )
  expect_false(grepl("np_beta_status", registry.header, fixed = TRUE))
  expect_false(grepl("np_beta_conditional_scalar_log(", native, fixed = TRUE))
  expect_false(grepl("NP_BETA_CONDITIONAL_", native, fixed = TRUE))
  expect_false(file.exists(file.path(root, "src", "beta_conditional.c")))
  expect_false(file.exists(file.path(root, "src", "beta_conditional.h")))
})
