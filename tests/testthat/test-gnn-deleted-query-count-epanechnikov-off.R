# Heavy numerical oracle: keep in a bounded fresh-pool singleton.
test_that("GNN deleted donor-count oracle: epanechnikov, trees off", {
  npRmpi_test_gnn_deleted_query_count(FALSE, "epanechnikov")
})
