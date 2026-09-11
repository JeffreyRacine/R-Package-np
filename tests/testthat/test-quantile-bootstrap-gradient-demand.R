test_that("quantile categorical demand retains the selected inverse pair", {
  pkg <- "npRmpi"
  helper <- getFromNamespace(".npqreg_categorical_first_differences", pkg)
  x <- data.frame(x = c(-.2, .1, .3),
                  u = factor(c("b", "a", "c")),
                  o = ordered(c("low", "middle", "high"),
                    levels = c("low", "middle", "high")))
  y <- data.frame(y = seq_len(nrow(x)))
  bws <- list(xndim = 3L, ixcon = c(TRUE, FALSE, FALSE),
              ixuno = c(FALSE, TRUE, FALSE), ixord = c(FALSE, FALSE, TRUE),
              xdati = list(icon = c(TRUE, FALSE, FALSE)))
  cache <- new.env(parent = emptyenv())
  cache$enabled <- FALSE
  calls <- list()
  local_mocked_bindings(
    .npqreg_invert_selected_cdf = function(bws, xdat, ydat, exdat, tau,
      tol, small, itmax, cdf.cache, cdf.row.keys, allow.external) {
      expect_identical(cdf.cache, cache)
      expect_false(allow.external)
      calls[[length(calls) + 1L]] <<- exdat
      exdat$x + 2 * as.integer(exdat$u) + 3 * as.integer(exdat$o)
    }, .package = pkg)
  args <- list(bws = bws, xdat = x, ydat = y, exdat = x,
    tau = .5, tol = .0001, small = .00001, itmax = 100L, cdf.cache = cache)
  full <- do.call(helper, args)
  full.calls <- calls
  expect_length(full.calls, 4L)
  for (j in 2:3) {
    calls <- list()
    one <- do.call(helper, c(args, list(gradient.target = j)))
    expect_identical(one[, j], full[, j])
    expect_identical(dim(one), dim(full))
    expect_identical(calls, full.calls[(2L * j - 3L):(2L * j - 2L)])
    expect_true(all(is.na(one[, -j, drop = FALSE])))
  }
  expect_error(do.call(helper, c(args, list(gradient.target = 2.5))),
               "invalid conditional gradient coordinate")
})

test_that("quantile scalar bootstrap demand preserves full-evaluator columns", {
  pkg <- "npRmpi"
  ns <- asNamespace(pkg)
  with.local <- if (pkg == "npRmpi") get(".npRmpi_with_local_regression", ns)
    else function(expr) force(expr)
  eval.raw <- get(".np_plot_quantile_eval", ns)
  eval.fun <- function(...) with.local(eval.raw(...))
  boot.raw <- get(".np_inid_boot_from_quantile_gradient_local", ns)
  boot.fun <- function(...) with.local(boot.raw(...))
  i <- seq_len(36L)
  x <- data.frame(x = sin(i * .37),
    u = factor(rep(c("a", "b"), length.out = length(i))),
    o = ordered(rep(c("low", "mid", "high"), length.out = length(i)),
      levels = c("low", "mid", "high")))
  y <- data.frame(y = sin(i * .31) + x$x)
  bws <- with.local(get("npcdistbw", ns)(xdat = x, ydat = y,
    bws = c(.8, .8, .2, .2), bandwidth.compute = FALSE, regtype = "ll"))
  ex <- x[c(11L, 27L), , drop = FALSE]
  counts <- rep.int(1L, nrow(x))
  counts[1:2] <- c(2L, 0L)
  idx <- rep.int(seq_len(nrow(x)), counts)
  args <- list(bws = bws, txdat = x, tydat = y, exdat = ex,
    tau = c(.4, .6), gradients = TRUE, need.errors = FALSE)
  full <- do.call(eval.fun, args)
  expanded <- args
  expanded$txdat <- x[idx, , drop = FALSE]
  expanded$tydat <- y[idx, , drop = FALSE]
  full.expanded <- do.call(eval.fun, expanded)
  for (j in seq_len(ncol(x))) {
    one <- do.call(eval.fun, c(args, list(gradient.target = j)))
    expect_identical(one$quantile, full$quantile)
    expect_identical(one$quantgrad[, j, ], full$quantgrad[, j, ])
    expect_identical(dim(one$quantgrad), dim(full$quantgrad))
    out <- boot.fun(xdat = x, ydat = y[[1L]], exdat = ex, bws = bws,
      B = 1L, tau = args$tau, gradient.index = j,
      counts = matrix(counts, ncol = 1L))
    expect_identical(out$t0, as.vector(full$quantgrad[, j, ]))
    expect_identical(as.vector(out$t), as.vector(full.expanded$quantgrad[, j, ]))
  }
  bad <- args
  bad$need.errors <- TRUE
  expect_error(do.call(eval.fun, c(bad, list(gradient.target = 1L))),
               "gradients=TRUE and need.errors=FALSE")
  bad <- args
  bad$gradients <- FALSE
  expect_error(do.call(eval.fun, c(bad, list(gradient.target = 1L))),
               "gradients=TRUE and need.errors=FALSE")
})
