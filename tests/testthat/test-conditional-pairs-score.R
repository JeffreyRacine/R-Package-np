test_that("conditional CDF pairs recalculate both denominators for resampled rows", {
  x <- c(-1, -.2, .6, 1.4)
  y <- c(.5, 1, -1, 2)
  grid <- data.frame(x = c(-.3, .5), y = c(.8, 1.6))
  h <- c(.4, .6)
  ids <- c(1L, 1L, 2L, 4L)
  ratio <- function(s, rows) {
    w <- dnorm(outer(x[rows], grid$x, "-") / (s*h[1]))
    num <- w * pnorm(-outer(y[rows], grid$y, "-") / (s*h[2]))
    list(w = w, num = num, value = colSums(num) / colSums(w))
  }
  a <- ratio(1, 1:4)
  b <- ratio(sqrt(2), 1:4)
  center <- 2*a$value - b$value
  plan <- list(ix = 1:2, w = a$w, jw = b$w, main = a$num, wide = b$num, center = center)
  got <- .np_conditional_pairs_block(list(value = matrix(tabulate(ids, 4), 4)), plan)
  expect_equal(as.numeric(got), 2*ratio(1, ids)$value - ratio(sqrt(2), ids)$value - center)
  plan$w[] <- 0
  expect_error(.np_conditional_pairs_block(list(value = matrix(1L, 4)), plan), "predictor mass")
})

test_that("conditional score and basic pairs companions retain signed endpoints and scopes", {
  control <- np_conditional_bias_control("pairs-score")
  expect_identical(control$construction,"pairs-score")
  g <- data.frame(x=c(-.2,.2),y=c(.1,.4))
  x <- seq(-1,1,length.out=31); y <- sin(x*3)
  B <- 199L
  roots <- rbind(seq(-2,4,length.out=B),seq(-.5,.9,length.out=B))
  local_mocked_bindings(.np_conditional_dispatch = function(tasks, plan, family, m) {
    expect_identical(family,"pairs-score"); roots
  }, .package = "npRmpi")
  run <- function() .np_with_seed(32L, .np_conditional_pairs_score(
    x, y, g, c(.4,.45), B, .05, rep(.5,2L), 0))
  z <- run()
  q <- t(apply(roots,1L,quantile,c(.025,.975),type=1,names=FALSE))
  expect_equal(z$bands$pointwise[1L,1L],unname(z$center[1L]-q[1L,2L]),ignore_attr=TRUE)
  expect_equal(z$bands$pointwise[1L,2L],unname(z$center[1L]-q[1L,1L]),ignore_attr=TRUE)
  expect_equal(z$critical.raw,unname(quantile(apply(abs(roots),2,max),.95,type=1)))
  expect_true(all(z$bands$simultaneous[,1]<=z$bands$pointwise[,1]))
  expect_true(all(z$bands$simultaneous[,2]>=z$bands$pointwise[,2]))
  roots[,] <- 0
  roots[1,1:4] <- 5
  rare <- run()
  expect_gt(rare$critical,rare$critical.raw)
  expect_equal(unname(rare$bands$bonferroni[1,1]),rare$center[1]-5)
  expect_identical(rare$law,"inid")
  label <- .np_plot_variability_descriptor("bootstrap","simultaneous",.05,"bias-corrected",construction="pairs-score")
  expect_true(any(grepl("score",unlist(label),fixed=TRUE)))
  expect_true(any(grepl("pairs",unlist(label),fixed=TRUE)))
})

test_that("signed CDF donor ranges account for interior extrema", {
  w <- c(.02,.1,.88); w2 <- c(.2,.4,.4)
  z <- list(w=matrix(w,3),jw=matrix(w2,3),ix=1L)
  b <- .np_conditional_pairs_range(z)
  u <- seq(-12,12,length.out=20001L)
  values <- outer(2*w,pnorm(u)) - outer(w2,pnorm(u/sqrt(2)))
  lo <- apply(values,1,min);hi <- apply(values,1,max)
  expect_lt(abs(b$lower-sum(lo)),2e-7)
  expect_lt(max(abs(b$span*b$q[,1]-(hi-lo))),2e-7)
})

test_that("the explicit conditional CDF plot records both interval components", {
  skip_if_not(.mpi_pool_active() && !.mpi_suite_local_mode_owned(), "requires the active distributed pool")
  x <- seq(-1.5, 1.5, length.out = 70)
  d <- data.frame(x=x, y=sin(3*x)+.3*cos(seq_along(x)))
  bw <- npcdistbw(y ~ x, data=d, bws=c(.33,.4), regtype="lc", bandwidth.compute=FALSE)
  fit <- npcdist(bw)
  ans <- plot(fit, center="bias-corrected", errors="bootstrap", bootstrap="inid",
    interval="percentile", bias.control=np_conditional_bias_control("pairs-score"),
    band="all", neval=3L, B=399L, output="data", random.seed=309L)
  info <- ans$inference[[1L]]
  expect_identical(info$resampling_law,"inid")
  expect_identical(info$correction,"pairs-score")
  expect_match(info$primary,"sharp-score")
  expect_match(info$endpoint_construction,"basic pairs")
  expect_gte(info$critical_value,info$critical_value_raw)
  expect_setequal(names(info$components),c("pointwise","simultaneous","bonferroni"))
  control <- np_conditional_bias_control("pairs-score")
  expect_warning(
    large <- plot(fit, center = "bias-corrected", errors = "bootstrap",
                  bootstrap = "inid", interval = "percentile", band = "all",
                  bias.control = control, neval = 26L, B = 399L,
                  output = "data", random.seed = 309L),
    "Bonferroni tails have limited finite-draw resolution")
  frame <- as.data.frame(large, panel = 1L, target = 1L, band = "pointwise")
  small <- as.data.frame(ans, panel = 1L, target = 1L, band = "pointwise")
  shared <- merge(small, frame, by = c("eval_1", "eval_2"), suffixes = c(".small", ".large"))
  expect_equal(nrow(frame), 26L * 26L)
  expect_gte(nrow(shared), 4L)
  expect_equal(shared$lower.small, shared$lower.large, tolerance = 1e-12)
  expect_equal(shared$upper.small, shared$upper.large, tolerance = 1e-12)
  expect_true(all(is.finite(frame$lower) & frame$lower <= frame$upper))
  expect_identical(large$inference[[1L]]$correction, "pairs-score")
  expect_error(plot(fit, center="bias-corrected", errors="bootstrap", bootstrap="geom",
    bias.control=np_conditional_bias_control("pairs-score"), band="simultaneous",
    neval=3L, B=399L, output="data"),"pairs-score|block")
})

test_that("only explicit pairs-score admits affordable larger conditional grids", {
  x <- seq(-1, 1, length.out = 80L)
  g <- expand.grid(x = x[1:17], y = x[1:17])
  expect_error(.np_conditional_workspace(x, x, g, c(.3, .4), 199L), "256 target")
  expect_lt(.np_conditional_workspace(x, x, g, c(.3, .4), 199L, TRUE)$estimated.bytes,
            256 * 1024^2)
  x.long <- seq(-1, 1, length.out = 1400L)
  for (ny in c(25L, 26L)) {
    rectangular <- expand.grid(x = x[1:10], y = x[seq_len(ny)])
    expect_lt(.np_conditional_workspace(x.long, x.long, rectangular,
                                        c(.3, .4), 399L, TRUE)$estimated.bytes,
              256 * 1024^2)
  }
  full <- expand.grid(x = x[1:50], y = x[1:50])
  expect_error(.np_conditional_workspace(x, x, full, c(.3, .4), 9999L, TRUE), "256 MiB")
  expect_error(.np_conditional_workspace(x, x, full[c(seq_len(2500), 1), ],
                                         c(.3, .4), 99L, TRUE), "2500 target")
  expect_error(.np_conditional_workspace(x, x, data.frame(x = 0, y = x[1:51]),
                                         c(.3, .4), 99L, TRUE), "50 distinct")
})
