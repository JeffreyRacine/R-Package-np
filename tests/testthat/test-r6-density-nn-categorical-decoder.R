test_that("NN density MADS evaluates decoded categorical coordinates", {
  skip_on_cran()
  if (!spawn_mpi_slaves(1L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  withr::local_options(np.messages = FALSE, np.extendednn = TRUE)
  withr::local_preserve_seed()
  set.seed(42)
  d <- data.frame(x = rnorm(40), z = factor(sample(letters[1:3],40,TRUE)))
  for (type in c("generalized_nn", "adaptive_nn")) {
    for (ordered in c(FALSE, TRUE)) {
      d$z <- factor(d$z, ordered = ordered)
      b <- npudensbw(dat = d, bwtype = type, bwsolver = "mads", nmulti = 1L,
                     nomad.opts = list(MAX_BB_EVAL = 50L))
      selected <- b$nomad.restart.results[[b$nomad.best.restart]]
      # Before the decoder port this fixture produces the 1e7 native penalty,
      # even though the endpoint's independent certificate looks plausible.
      expect_lt(abs(selected$native$objective), 1e6)
      expect_equal(selected$native$objective, selected$objective,
                   tolerance = 2e-12)
      expect_identical(as.numeric(b$fval), -as.numeric(selected$objective))
    }
  }
})
