npRmpi_full_test_singleton_files <- function() {
  c(
    "test-regression-hc0-general-lp-derivative-contract.R",
    "test-regression-hc0-scalar-categorical-contract.R",
    "test-session-routing-subprocess-contract.R"
  )
}

npRmpi_full_test_plan <- function(test_files, shard_size) {
  if (!is.character(test_files) || !length(test_files) ||
      anyNA(test_files) || any(!nzchar(test_files)) ||
      anyDuplicated(test_files)) {
    stop("npRmpi full-suite plan requires unique test files")
  }
  if (!identical(test_files, sort(test_files)))
    stop("npRmpi full-suite plan requires sorted test files")
  shard_size <- suppressWarnings(as.integer(shard_size))
  if (length(shard_size) != 1L || is.na(shard_size) || shard_size < 1L)
    stop("npRmpi full-suite plan requires a positive shard size")

  singleton_files <- npRmpi_full_test_singleton_files()
  if (!is.character(singleton_files) || anyNA(singleton_files) ||
      any(!nzchar(singleton_files)) || anyDuplicated(singleton_files) ||
      !all(singleton_files %in% test_files)) {
    stop("npRmpi full-suite singleton registry is invalid")
  }

  plan <- list()
  pending <- character()
  publish <- function(files) {
    plan[[length(plan) + 1L]] <<- files
  }

  for (test_file in test_files) {
    if (test_file %in% singleton_files) {
      if (length(pending)) {
        publish(pending)
        pending <- character()
      }
      publish(test_file)
    } else {
      pending <- c(pending, test_file)
      if (length(pending) == shard_size) {
        publish(pending)
        pending <- character()
      }
    }
  }
  if (length(pending))
    publish(pending)

  flattened <- unlist(plan, use.names = FALSE)
  if (!identical(flattened, test_files) ||
      any(lengths(plan) < 1L) || any(lengths(plan) > shard_size) ||
      !all(vapply(singleton_files, function(test_file) {
        any(vapply(plan, identical, logical(1), test_file))
      }, logical(1)))) {
    stop("npRmpi full-suite plan lost or duplicated a test file")
  }
  plan
}
