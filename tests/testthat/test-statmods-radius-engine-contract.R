library(npRmpi)

locate_statmods_c <- function() {
  candidates <- c(
    test_path("..", "..", "src", "statmods.c"),
    test_path("..", "..", "..", "src", "statmods.c"),
    file.path(Sys.getenv("R_PACKAGE_DIR", ""), "src", "statmods.c"),
    file.path(Sys.getenv("R_PACKAGE_SOURCE", ""), "src", "statmods.c"),
    file.path(getwd(), "src", "statmods.c"),
    file.path(getwd(), "..", "src", "statmods.c")
  )
  candidates <- unique(candidates[nzchar(candidates)])
  hits <- candidates[file.exists(candidates)]
  if (length(hits) == 0L) {
    return(NULL)
  }
  hits[[1L]]
}

test_that("statmods uses observation-count support for adaptive and generalized radii", {
  src_file <- locate_statmods_c()
  skip_if(is.null(src_file), "source file src/statmods.c unavailable in this test context")

  lines <- readLines(src_file, warn = FALSE)

  expect_true(any(grepl("compute_nn_distance_observation_support_subset", lines, fixed = TRUE)))
  expect_true(any(grepl("compute_nn_distance_train_eval_observation_support_subset", lines, fixed = TRUE)))
  expect_true(any(grepl("kth_observation_radius_from_support", lines, fixed = TRUE)))
  expect_true(any(grepl("kth_observation_radius_for_eval_from_support", lines, fixed = TRUE)))
  expect_true(any(grepl("nearest_positive_radius_from_support", lines, fixed = TRUE)))

  expect_equal(sum(grepl("sort_safe(num_obs_train, vector_dist)", lines, fixed = TRUE)), 0L)
  expect_equal(sum(grepl("vector_unique_dist", lines, fixed = TRUE)), 0L)
  expect_equal(sum(grepl("pointer_nndi", lines, fixed = TRUE)), 0L)

  expect_gte(sum(grepl("compute_nn_distance_observation_support_subset\\(", lines)), 3L)
  expect_gte(sum(grepl("compute_nn_distance_train_eval_observation_support_subset\\(", lines)), 3L)
})
