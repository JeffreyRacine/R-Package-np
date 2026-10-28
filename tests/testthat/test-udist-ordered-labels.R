test_that("GNN corrected centers retain canonical ordered label distances", {
  x <- data.frame(x = sin(1:40), y = cos(1:40),
    o = ordered(rep(c("low", "middle", "high", "top"), 10),
      levels = c("low", "middle", "high", "top", "unused")))
  ranked <- x; ranked$o <- ordered(as.integer(x$o), levels = 1:5)
  make <- function(data) {
    bw <- npudistbw(dat = data, bws = c(12, 13, .2),
      bwtype = "generalized_nn", bandwidth.compute = FALSE)
    .np_udist_bc_prepare(bw, data)
  }
  state <- make(x); reference <- make(ranked)
  grid <- x[c(1, 9, 22, 1, 40), ]; rg <- ranked[c(1, 9, 22, 1, 40), ]
  expect_identical(.np_udist_gnn_main(state, grid),
    .np_udist_gnn_main(reference, rg))
  a <- .np_udist_gnn_blocks(state, grid, 2)
  b <- .np_udist_gnn_blocks(reference, rg, 2)
  for (j in seq_along(a)) {
    pa <- .np_udist_gnn_block(state, a[[j]])
    pb <- .np_udist_gnn_block(reference, b[[j]])
    expect_identical(pa, pb)
  }
})
