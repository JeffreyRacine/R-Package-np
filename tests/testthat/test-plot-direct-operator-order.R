test_that("direct plot weights keep operators with their data columns", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  dat <- data.frame(f = factor(rep(c("a", "b", "c"), 8L)),
                    x = seq(-1, 1, length.out = 24L),
                    o = ordered(rep(c("l", "h"), 12L)),
                    z = sin(seq_len(24L)))
  bandwidths <- c(f = .3, x = .5, o = .4, z = .6)
  for (columns in list(c("f", "x", "o", "z"), c("z", "o", "x", "f"),
                       c("x", "z", "f", "o"))) {
    td <- dat[, columns]
    ed <- td[c(3L, 10L, 20L), , drop = FALSE]
    bw <- npudensbw(dat = td, bws = bandwidths[columns], bandwidth.compute = FALSE)
    for (operator in list(rep("normal", 4L), rep("integral", 4L),
                          c(f = "normal", x = "normal", o = "normal", z = "integral")[columns],
                          c(f = "normal", x = "derivative", o = "normal", z = "normal")[columns])) {
      got <- npRmpi:::.np_plot_kernel_weights_direct(bw, td, ed, operator)
      expected <- npRmpi:::.np_kernel_weights_direct(bw, td, ed, operator = operator)
      expect_equal(got, expected, tolerance = 1e-14)
    }
  }
})
