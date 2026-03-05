library(npRmpi)

test_that("CVLS LL/LP route predicate is centralized in one helper", {
  src_file <- test_path("..", "..", "src", "jksum.c")
  expect_true(file.exists(src_file))

  lines <- readLines(src_file, warn = FALSE)

  expect_true(any(grepl("np_reg_cv_use_symmetric_dropone_path", lines, fixed = TRUE)))

  raw_pred <- "(bwm == RBWM_CVLS) || ks_tree_use || (BANDWIDTH_reg == BW_ADAP_NN)"
  # Exactly one raw predicate instance is allowed (the helper definition itself).
  expect_equal(sum(grepl(raw_pred, lines, fixed = TRUE)), 1L)
  # The two-term suffix should also only appear inside the helper definition.
  expect_equal(sum(grepl("ks_tree_use || (BANDWIDTH_reg == BW_ADAP_NN)", lines, fixed = TRUE)), 1L)
  # No inline ad hoc route predicates should remain in if-statements.
  expect_equal(sum(grepl("if\\(\\(bwm == RBWM_CVLS\\)", lines)), 0L)
  expect_equal(sum(grepl("if\\(ks_tree_use \\|\\| \\(BANDWIDTH_reg == BW_ADAP_NN\\)\\)", lines)), 0L)

  helper_calls <- sum(grepl("np_reg_cv_use_symmetric_dropone_path\\(", lines))
  expect_gte(helper_calls, 3L)
})
