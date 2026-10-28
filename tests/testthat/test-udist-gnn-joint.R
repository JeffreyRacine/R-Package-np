test_that("exact GNN centers preserve expanded samples and requested target pairs", {
  x <- data.frame(x = sin(1:31), y = cos(1:31),
    o = ordered(rep(c(-2, 0, 3), length.out = 31)))
  bw <- npudistbw(dat = x, bws = c(10, 11, .2),
    bwtype = "generalized_nn", bandwidth.compute = FALSE)
  state <- .np_udist_bc_prepare(bw, x)
  grid <- x[c(7, 2, 11, 2, 1, 31), , drop = FALSE]
  set.seed(610271L)
  counts <- as.vector(rmultinom(1L, 31L, rep(1/31, 31)))
  index <- rep.int(1:31, counts)
  sample <- state; sample$x <- x[index, , drop = FALSE]
  sample$category[[1L]]$training <- state$category[[1L]]$training[index]
  reference <- .np_udist_gnn_main(sample, grid)$center
  answer <- numeric(nrow(grid)); covered <- integer()
  for (block in .np_udist_gnn_blocks(state, grid, size = 2L)) {
    prepared <- .np_udist_gnn_block(state, block)
    answer[block$rows] <- .np_udist_gnn_center(prepared, counts)$center
    covered <- c(covered, block$rows)
  }
  expect_identical(sort(covered), seq_len(nrow(grid)))
  expect_equal(answer, reference, tolerance = 1e-12)
  expect_identical(answer[2L], answer[4L])
})

test_that("GNN absolute bands retain nesting, singleton equality and RNG", {
  old <- options(np.messages = FALSE, np.plot.progress = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = sin(1:31))
  bw <- npudistbw(dat = x, bws = 10, bwtype = "generalized_nn",
    bandwidth.compute = FALSE)
  state <- .np_udist_bc_prepare(bw, x)
  for (law in c("gaussian", "inid", "fixed", "geom")) {
    set.seed(610271L)
    result <- .np_udist_joint_inference(state, x[1:3, , drop = FALSE],
      law = law, B = 49L, alpha = c(.05, .1),
      blocklen = if (law %in% c("fixed", "geom")) 3L else NULL)
    iv <- result$intervals
    expect_true(all(iv$pointwise[, 1L, ] >= iv$simultaneous[, 1L, ]))
    expect_true(all(iv$pointwise[, 2L, ] <= iv$simultaneous[, 2L, ]))
    expect_true(all(iv$simultaneous[, 1L, ] >= iv$bonferroni[, 1L, ]))
    expect_true(all(iv$simultaneous[, 2L, ] <= iv$bonferroni[, 2L, ]))
    expect_identical(result$rng$end, .Random.seed)
    one <- .np_udist_joint_inference(state, x[1L, , drop = FALSE],
      law = law, B = 49L, blocklen = if (law %in% c("fixed", "geom")) 3L else NULL)
    expect_identical(one$intervals$pointwise, one$intervals$simultaneous)
    expect_identical(one$intervals$simultaneous, one$intervals$bonferroni)
  }
})

test_that("zero GNN radii are unavailable and memory partitions are explicit", {
  x <- data.frame(x = seq_len(31))
  bw <- npudistbw(dat = x, bws = 5, bwtype = "generalized_nn",
    bandwidth.compute = FALSE)
  state <- .np_udist_bc_prepare(bw, x)
  block <- .np_udist_gnn_blocks(state, x[1L, , drop = FALSE])[[1L]]
  expect_error(.np_udist_gnn_center(.np_udist_gnn_block(state, block),
    c(31L, rep(0L, 30L))), "nonpositive physical radius")
  state$max.bytes <- 1
  expect_error(.np_udist_gnn_budget(state, 1L, 49L), "memory forecast")
})
