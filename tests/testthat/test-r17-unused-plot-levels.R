test_that("categorical bootstrap summaries follow evaluated labels, not observed levels", {
  summarise <- getFromNamespace(".np_plot_boot_factor_boxplots", "npRmpi")
  describe <- getFromNamespace("untangle", "npRmpi")
  f <- factor(c("b", "d", "b"), levels = c("a", "b", "c", "d", "e"))
  td <- describe(data.frame(f = f))
  draws <- matrix(seq_len(21), 7L, 3L)
  targets <- factor(c("e", "b", "a"), levels = levels(f))
  bp <- summarise(draws, td, 1L, 7L, eval.values = targets)
  expect_identical(bp$names, c("e", "b", "a"))
  expect_equal(bp$stats, vapply(seq_len(3L), function(j) boxplot.stats(draws[, j])$stats, numeric(5L)))
  expect_error(summarise(draws, td, 1L, 7L, eval.values = targets[1L]), "evaluation")
  expect_error(summarise(draws, td, 1L, 7L, eval.values = c("e", "b", "z")), "evaluation")
})

test_that("conditional and semiparametric plots share categorical draw labels", {
  set.seed(18118)
  n <- 48L
  d <- data.frame(x = runif(n), z = runif(n),
                  u = factor(rep(c("b", "c"), n / 2L), levels = c("a", "b", "c")))
  d$y <- d$x + as.integer(d$u) + rnorm(n, sd = .2)
  models <- list(
    cdens = npcdens(y ~ x + u, data = d, bws = c(.5, .5, .2), bandwidth.compute = FALSE),
    cdist = npcdist(y ~ x + u, data = d, bws = c(.5, .5, .2), bandwidth.compute = FALSE),
    plreg = npplreg(y ~ x | z + u, data = d, bws = matrix(c(.5, .5, .2, .2), 2L),
                   bandwidth.compute = FALSE),
    scoef = npscoef(y ~ x | u, data = d, bws = .2, bandwidth.compute = FALSE,
                   iterate = FALSE))
  boxes <- function(x) {
    if (!is.list(x)) return(list())
    if (all(c("stats", "conf", "n", "names") %in% names(x))) return(list(x))
    unlist(lapply(x, boxes), recursive = FALSE)
  }
  for (name in names(models)) {
    for (gradient in if (name %in% c("cdens", "cdist")) c(FALSE, TRUE) else FALSE) {
      args <- list(x = models[[name]], errors = "bootstrap", B = 7L, neval = 3L,
                   band = "all", perspective = FALSE, common.scale = FALSE,
                   plot.behavior = "data")
      if (name != "plreg") args$gradients <- gradient
      out <- suppressWarnings(do.call(plot, args))
      bp <- boxes(out)
      expect_true(length(bp) > 0L, info = name)
      for (box in bp) {
        expect_identical(box$names, levels(d$u), info = name)
        expect_equal(ncol(box$stats), 3L, info = name)
      }
    }
  }
})

test_that("unused factor levels remain evaluated and labelled in bootstrap plots", {
  set.seed(17118)
  n <- 48L
  x <- runif(n)
  y <- x + rep(c(0, 1), n / 2L) + rnorm(n, sd = .2)
  for (ordered in c(FALSE, TRUE)) for (position in seq_len(3L)) {
    # Numeric ordered labels describe lattice coordinates, not a permutation
    # of ranks. Move the unused observation level, not that ordered lattice.
    lev <- if (ordered) c("0", "1", "2") else
      list(c("0", "1", "2"), c("1", "0", "2"), c("1", "2", "0"))[[position]]
    observed <- if (ordered) lev[-position] else c("1", "2")
    f <- factor(rep(observed, n / 2L), levels = lev, ordered = ordered)
    d <- data.frame(y = y, x = x, f = f)
    for (fam in c("reg", "dens", if (ordered) "dist")) {
      b <- switch(fam,
        reg = npregbw(y ~ x + f, data = d, bws = c(.4, .2), bandwidth.compute = FALSE),
        dens = npudensbw(~ x + f, data = d, bws = c(.4, .2), bandwidth.compute = FALSE),
        dist = npudistbw(~ x + f, data = d, bws = c(.4, .2), bandwidth.compute = FALSE))
      set.seed(7118)
      out <- suppressWarnings(plot(b, plot.behavior = "data", errors = "bootstrap", B = 7,
                  neval = 5L, common.scale = FALSE, band = "all", perspective = FALSE))
      expect_identical(out[[2L]]$bxp$names, lev)
      expect_equal(ncol(out[[2L]]$bxp$stats), 3L)
      expect_true(all(is.finite(out[[2L]]$bxp$stats)))
      expect_identical(levels(d$f), lev)
    }
  }
})
