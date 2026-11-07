# Unique default-lane contract: refuse non-representable normal task state
# before advance. No pool, native fixture or heavy bootstrap is required.
test_that("seed-only normal tasks refuse Box-Muller before advancement", {
  kind <- RNGkind(); withr::defer(do.call(RNGkind, as.list(kind)))
  withr::local_preserve_seed()
  RNGkind("Mersenne-Twister", "Box-Muller", "Rejection"); set.seed(511)
  attach <- get(".npRmpi_bootstrap_attach_rng_stream", asNamespace("npRmpi"))
  for (pending in c(FALSE, TRUE)) {
    set.seed(511); if (pending) rnorm(1)
    before <- .Random.seed; calls <- 0L
    for (n in c(1L, 2L)) {
      tasks <- rep(list(list(start=1L,bsz=1L)), n)
      expect_error(attach(tasks, function(task) { calls <<- calls+1L; rnorm(1) },
                          normal.draws=TRUE), "seed-only bootstrap tasks cannot reproduce Box-Muller", fixed=TRUE)
      expect_identical(.Random.seed,before); expect_identical(calls,0L)
    }
  }
  expect_identical(attach(list(), function(task) stop("advanced"), TRUE),list())
  expect_error(attach(list(list()), NULL, TRUE), "advance must be a function")
})

test_that("normal-free task streams retain exact Box-Muller uniforms", {
  kind <- RNGkind(); withr::defer(do.call(RNGkind, as.list(kind)))
  withr::local_preserve_seed()
  RNGkind("Mersenne-Twister", "Box-Muller", "Rejection"); set.seed(533)
  attach <- get(".npRmpi_bootstrap_attach_rng_stream", asNamespace("npRmpi"))
  tasks <- list(list(bsz=3L),list(bsz=2L)); before <- .Random.seed
  wanted <- runif(5); end <- .Random.seed; assign(".Random.seed",before,.GlobalEnv)
  values <- numeric(); answer <- attach(tasks,function(task) {values <<- c(values,runif(task$bsz))})
  expect_identical(values,wanted);expect_identical(.Random.seed,end)
  expect_identical(attr(answer,"rng_final_state"),end)
})
