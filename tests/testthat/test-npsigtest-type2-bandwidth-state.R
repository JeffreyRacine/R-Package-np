test_that("type-II replacement preserves every bandwidth representation", {
  replace <- getFromNamespace(".np_npsig_replace_bandwidth", "np")
  physical <- getFromNamespace(".np_physical_bandwidth", "np")
  set.seed(169)
  n <- 32L
  x <- data.frame(u = factor(rep(c("a", "b"), n/2)), x = runif(n),
                  o = ordered(rep(1:4, n/4)), z = rnorm(n))
  y <- x$x + rnorm(n)
  for (type in c("fixed", "generalized_nn", "adaptive_nn"))
    for (scaled in c(FALSE, TRUE))
      for (new.scaled in c(FALSE, TRUE)) {
        old.bw <- c(.2, if (type == "fixed") .7 else 12, .3,
                    if (type == "fixed") .8 else 14)
        new.bw <- c(.1, if (type == "fixed") .4 else 13, .2,
                    if (type == "fixed") 1.1 else 15)
        original <- npregbw(xdat = x, ydat = y, bws = old.bw,
          bwtype = type, bwscaling = scaled, bandwidth.compute = FALSE)
        selected <- npregbw(xdat = transform(x, x = x * 2, z = z * 3),
          ydat = y, bws = new.bw, bwtype = type, bwscaling = new.scaled,
          bandwidth.compute = FALSE)
        saved <- serialize(list(original, selected), NULL)
        seed <- .Random.seed
        for (index in list(1L, 2L, 3L, 4L, c(1L, 4L), 1:4)) {
          out <- replace(original, selected, index)
          expected <- physical(original)
          expected[index] <- physical(selected)[index]
          expect_identical(unname(physical(out)), unname(expected))
          for (field in c("bandwidth", "sfactor", "sumNum")) {
            target <- original[[field]][["x"]]
            target[index] <- selected[[field]][["x"]][index]
            expect_identical(out[[field]][["x"]], target)
          }
          expect_identical(out$scaling, original$scaling)
          expect_identical(out$bw[-index], original$bw[-index])
          target.sd <- original$sdev
          positions <- which(which(original$icon) %in% index)
          target.sd[positions] <- selected$sdev[positions]
          expect_identical(out$sdev, target.sd)
        }
        expect_identical(serialize(list(original, selected), NULL), saved)
        expect_identical(.Random.seed, seed)
      }
})

test_that("type-II scaled tests match physical oracles and leave type I alone", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(1691)
  x <- data.frame(x = runif(32), z = rnorm(32))
  y <- sin(3*x$x) + x$z + rnorm(32, sd = .3)
  original <- npregbw(xdat = x, ydat = y, bws = c(.9, 1.3),
    bwscaling = TRUE, bandwidth.compute = FALSE)
  physical <- npregbw(xdat = x, ydat = y, bws = original$bandwidth$x,
    bandwidth.compute = FALSE)
  selected <- npregbw(xdat = x, ydat = y, bws = c(1.4, .6),
    bwscaling = TRUE, bandwidth.compute = FALSE)
  selected.physical <- npregbw(xdat = x, ydat = y,
    bws = selected$bandwidth$x, bandwidth.compute = FALSE)
  state <- new.env(parent = emptyenv())
  state$selected <- selected
  state$count <- 0L
  stub <- function(...) {
    state$count <- state$count + 1L
    state$selected
  }
  local_mocked_bindings(.np_npsig_bootstrap_bw_reselect = stub, .package = "np")
  saved <- serialize(list(original, physical, selected, selected.physical), NULL)
  seed <- .Random.seed
  cells <- list(list(index = 2L, joint = FALSE),
                list(index = 1:2, joint = FALSE),
                list(index = 1:2, joint = TRUE))
  for (boot.method in c("iid", "wild", "wild-rademacher", "pairwise"))
    for (cell in cells) {
      state$selected <- selected
      state$count <- 0L
      a <- do.call(npsigtest, c(list(bws = original, xdat = x, ydat = y,
                    B = 9L, boot.type = "II", boot.method = boot.method), cell))
      expected <- 9L * if (cell$joint) 1L else length(cell$index)
      expect_identical(state$count, expected)
      state$selected <- selected.physical
      state$count <- 0L
      b <- do.call(npsigtest, c(list(bws = physical, xdat = x, ydat = y,
                    B = 9L, boot.type = "II", boot.method = boot.method), cell))
      expect_identical(state$count, expected)
      for (field in c("In", "In.bootstrap", "P"))
        expect_equal(a[[field]], b[[field]], tolerance = 1e-12)
    }
  state$count <- 0L
  npsigtest(original, xdat = x, ydat = y, B = 9L, boot.type = "I")
  expect_identical(state$count, 0L)
  expect_identical(serialize(list(original, physical, selected, selected.physical), NULL), saved)
  expect_identical(.Random.seed, seed)
})
