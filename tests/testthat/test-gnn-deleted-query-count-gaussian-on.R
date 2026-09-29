# Heavy numerical oracle: keep in a bounded fresh-pool singleton.
test_that("GNN deleted donor-count oracle: gaussian, trees on", {
  npRmpi_test_gnn_deleted_query_count(TRUE, "gaussian")
})
