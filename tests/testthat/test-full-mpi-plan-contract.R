test_that("the bounded MPI plan preserves every sorted test and heavy singleton", {
  helper <- testthat::test_path("..", "validation", "full_mpi_test_plan.R")
  skip_if_not(file.exists(helper), "source plan helper unavailable")
  scope <- new.env(parent = globalenv())
  sys.source(helper, envir = scope)
  files <- sort(list.files(testthat::test_path(), "^test-.*\\.[rR]$"))
  singletons <- scope$npRmpi_full_test_singleton_files()
  for (size in c(1L, 3L, 10L)) {
    plan <- scope$npRmpi_full_test_plan(files, size)
    expect_identical(unlist(plan, use.names = FALSE), files)
    expect_true(all(lengths(plan) >= 1L & lengths(plan) <= size))
    for (file in singletons)
      expect_equal(sum(vapply(plan, identical, TRUE, file)), 1L, info = file)
  }
  expect_error(scope$npRmpi_full_test_plan(rev(files), 10L), "sorted test files")
  expect_error(scope$npRmpi_full_test_plan(c(files, files[1L]), 10L), "unique test files")
  expect_error(scope$npRmpi_full_test_plan(files, 0L), "positive shard size")
  expect_error(scope$npRmpi_full_test_plan(setdiff(files, singletons[1L]), 10L),
               "singleton registry is invalid")
})

test_that("formula training ownership has an independent MPI budget", {
  helper <- testthat::test_path("..", "validation", "full_mpi_test_plan.R")
  skip_if_not(file.exists(helper), "source plan helper unavailable")
  scope <- new.env(parent = globalenv())
  sys.source(helper, envir = scope)
  required <- "test-formula-training-ownership.R"
  files <- sort(list.files(testthat::test_path(), "^test-.*\\.[rR]$"))
  expect_true(required %in% files)
  expect_true(required %in% scope$npRmpi_full_test_singleton_files())
  plan <- scope$npRmpi_full_test_plan(files, 10L)
  expect_equal(sum(vapply(plan, identical, TRUE, required)), 1L)
})

test_that("expensive conditional GNN oracles have independent MPI budgets", {
  helper <- testthat::test_path("..", "validation", "full_mpi_test_plan.R")
  skip_if_not(file.exists(helper), "source plan helper unavailable")
  scope <- new.env(parent = globalenv())
  sys.source(helper, envir = scope)
  # Deliberately independent of the registry: removing an entry must fail.
  required <- c("test-conditional-gnn-general.R",
                "test-conditional-gnn-prefix.R",
                "test-conditional-gnn-projected.R")
  files <- sort(list.files(testthat::test_path(), "^test-.*\\.[rR]$"))
  expect_true(all(required %in% files))
  expect_true(all(required %in% scope$npRmpi_full_test_singleton_files()))
  # Include ordinary neighbors to test both pending-group flush boundaries.
  synthetic <- sort(unique(c(scope$npRmpi_full_test_singleton_files(), required,
    "test-conditional-gnn-general-before.R", "test-conditional-gnn-general0.R",
    "test-conditional-gnn-prefix-before.R", "test-conditional-gnn-prefix0.R",
    "test-conditional-gnn-projected-before.R", "test-conditional-gnn-projected0.R")))
  for (inventory in list(files, synthetic)) {
    for (size in c(1L, 3L, 10L, length(inventory) + 1L)) {
      plan <- scope$npRmpi_full_test_plan(inventory, size)
      expect_identical(unlist(plan, use.names = FALSE), inventory)
      expect_true(all(lengths(plan) >= 1L & lengths(plan) <= size))
      for (file in required)
        expect_equal(sum(vapply(plan, identical, TRUE, file)), 1L, info = file)
    }
  }
})

test_that("the GNN singleton wrappers cover every explicit tree and kernel", {
  expected <- expand.grid(kernel = c("gaussian", "epanechnikov", "beta"),
                          mode = c("off", "on"), stringsAsFactors = FALSE)
  actual <- sort(list.files(testthat::test_path(),
    "^test-gnn-deleted-query-count-(gaussian|epanechnikov|beta)-(off|on)\\.R$"))
  names <- paste0("test-gnn-deleted-query-count-", expected$kernel, "-", expected$mode, ".R")
  expect_identical(actual, sort(names))
  for (i in seq_len(nrow(expected))) {
    expr <- parse(testthat::test_path(names[i]), keep.source = FALSE)
    expect_length(expr, 1L)
    call <- expr[[1L]][[3L]][[2L]]
    expect_identical(call, as.call(list(as.name("npRmpi_test_gnn_deleted_query_count"),
      expected$mode[i] == "on", expected$kernel[i])))
  }
  body <- body(npRmpi_test_gnn_deleted_query_count)
  loops <- Filter(function(x) is.call(x) && identical(x[[1L]], as.name("for")), as.list(body))
  expect_length(loops, 1L)
  expect_identical(loops[[1L]][[3L]], as.name("trees"))
  expect_identical(loops[[1L]][[4L]][[3L]], as.name("kernels"))
})
