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
  }, .package = "np")
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
  allocation <- .np_conditional_workspace(x, x, full, c(.3, .4), 9999L, TRUE)
  expect_lt(allocation$tile.size, nrow(full))
  expect_lte(allocation$estimated.bytes, 256 * 1024^2)
  expect_error(.np_conditional_workspace(x.long, x.long, full, c(.3, .4), 9999L, TRUE),
               "256 MiB")
  expect_error(.np_conditional_workspace(x, x, full[c(seq_len(2500), 1), ],
                                         c(.3, .4), 99L, TRUE), "2500 target")
  expect_error(.np_conditional_workspace(x, x, data.frame(x = 0, y = x[1:51]),
                                         c(.3, .4), 99L, TRUE), "50 distinct")
})

test_that("pairs-score target tiles preserve the complete finite-grid calculation", {
  x <- seq(-1.5, 1.5, length.out = 30L)
  y <- sin(x) + cos(seq_along(x)) / 3
  grid <- expand.grid(x = c(.3, -.6, .1), y = c(.4, -.2))
  grid <- grid[c(4L, 1L, 3L, 4L, 2L, 6L), ]
  run <- function(tile) .np_with_seed(309L, .np_conditional_pairs_score(
    x, y, grid, c(.4, .3), 129L, .05, rep(0, nrow(grid)), 0, tile))
  baseline <- run(NULL)
  for (size in c(1L, 4L)) {
    tiled <- run(size)
    for (field in c("center", "bands", "se", "critical.raw", "critical"))
      expect_equal(tiled[[field]], baseline[[field]], tolerance = 1e-12)
  }
})

test_that("the over-budget pairs selector limits both storage and transmissions", {
  choose <- function(n, nx, ny, draws) {
    x <- seq(-1, 1, length.out = n)
    g <- expand.grid(x = seq_len(nx), y = seq_len(ny))
    .np_conditional_workspace(x, x, g, c(.3, .4), draws, TRUE)
  }
  expect_identical(choose(400L, 50L, 50L, 9999L)$tile.size, 250L)
  expect_identical(choose(1000L, 18L, 18L, 9999L)$tile.size, 139L)
  expect_identical(choose(1400L, 18L, 18L, 9999L)$tile.size, 31L)
  expect_error(choose(1450L, 18L, 18L, 9999L), "budget")
  expect_error(choose(150L, 50L, 50L, 10000L), "9,999")
  expect_error(choose(557L, 50L, 50L, 9999L), "budget")
})

test_that("zero predictor mass in a late tile refuses before dispatch", {
  x <- c((1:9)/100, 20)
  y <- seq_along(x)
  grid <- data.frame(x = c(.03, .06, 20), y = c(2, 5, 1))
  local_mocked_bindings(.np_conditional_dispatch = function(...)
    stop("unexpected dispatch before checking support"), .package = "np")
  run <- function(g, tile = 1L) .np_conditional_pairs_score(
    x, y, g, c(.01, .4), 129L, .05, numeric(nrow(g)), 0, tile)
  for (tile in list(0L, 2.5, NA_real_, Inf, 4L))
    expect_error(run(grid, tile), "tile size")
  set.seed(123)
  state <- .Random.seed
  g <- grid
  g$x[3L] <- 1000
  expect_error(run(g), "undefined predictor mass")
  expect_identical(.Random.seed, state)
  expect_error(.np_with_seed(309L, run(grid)), "resampled conditional CDF has zero")
  expect_identical(.Random.seed, state)
})

test_that("the Epanechnikov response uses its native scale and exact signed range", {
  s <- sqrt(5)
  u <- c(-1.1, -1, -1 + 1e-8, -.25, 0, .25, 1 - 1e-8, 1, 1.1)
  reference <- vapply(u, function(t) {
    if (t <= -1) return(0)
    integrate(function(v) .75 * (1 - v^2), -1, min(t, 1))$value
  }, numeric(1))
  expect_equal(.np_conditional_epan_cdf(s * u), reference, tolerance = 1e-14)
  expect_gt(.np_conditional_epan_cdf(s * u[3]), 0)
  w <- c(0, .03, .6, .37, 0)
  jw <- c(.07, .5, .2, 0, .23)
  z <- list(w = matrix(w), jw = matrix(jw), ix = 1L)
  scores <- .np_conditional_pairs_range(z, "epanechnikov")
  u <- seq(-s * sqrt(2), s * sqrt(2), length.out = 30001L)
  primitive <- function(v) {
    v <- pmax(-1, pmin(1, v/s))
    .5 + .25*v*(3 - v^2)
  }
  signed <- outer(2*w, primitive(u)) - outer(jw, primitive(u / sqrt(2)))
  low <- apply(signed, 1L, min)
  high <- apply(signed, 1L, max)
  expect_equal(scores$lower, sum(low), tolerance = 1e-7)
  expect_equal(drop(scores$span*scores$q), high-low, tolerance = 1e-7)
})

test_that("the pairs companion admits only the requested compact response", {
  d <- data.frame(x = seq(-1, 1, length.out = 55L))
  d$y <- sin(2*d$x) + cos(seq_len(nrow(d)))
  fit <- function(kx, ky) {
    bw <- npcdistbw(y ~ x, data = d, bws = c(.35, .4), regtype = "lc",
                    bandwidth.compute = FALSE, cxkertype = kx, cykertype = ky)
    npcdist(bw)
  }
  a <- list(x = fit("gaussian", "epanechnikov"), center = "bias-corrected",
            errors = "bootstrap", bootstrap = "inid", interval = "percentile",
            band = "all", bias.control = np_conditional_bias_control("pairs-score"),
            neval = 3L, B = 399L, output = "data", random.seed = 309L)
  old <- RNGkind()
  # The public caller supplies a scoped seed and must restore its entry state.
  set.seed(45)
  state <- .Random.seed
  z <- do.call(plot, a)
  expect_identical(.Random.seed, state)
  expect_identical(RNGkind(), old)
  expect_identical(z$inference[[1L]]$correction, "pairs-score")
  v <- as.data.frame(z, panel = 1L, target = 1L, band = "simultaneous")
  uncorrected <- as.data.frame(plot(a$x, errors = "none", neval = 3L,
                                    output = "data"), panel = 1L, target = 1L)
  expect_equal(v$estimate, uncorrected$estimate, tolerance = 2e-8)
  expect_true(all(is.finite(v$lower) & v$lower <= v$upper))
  a$x <- fit("epanechnikov", "gaussian")
  expect_error(do.call(plot, a), "Gaussian2 X")
})
test_that("compact response pairs preserve endpoint order, saturation and ties", {
  x <- seq(-1, 1, length.out = 48L)
  y <- rep(c(-.8, 0, .8), each = 16L)
  g <- data.frame(x = c(0, -.2, .2, -.2, 0), y = c(9, 0, -9, 0, 9))
  ord <- c(3L, 5L, 1L, 4L, 2L)
  draw <- function(g) .np_with_seed(312L, .np_conditional_pairs_score(
    x, y, g, c(.43, .18), 199L, .05, rep(0, nrow(g)), 0,
    response.kernel = "epanechnikov"))
  z <- draw(g)
  zz <- draw(g[ord, , drop = FALSE])
  expect_equal(z$center[c(1L, 3L, 5L)], c(1, 0, 1), tolerance = 1e-14)
  expect_equal(z$center[ord], zz$center, tolerance = 1e-14)
  for (sc in names(z$bands)) {
    expect_equal(z$bands[[sc]][ord, , drop = FALSE], zz$bands[[sc]], tolerance = 1e-14)
    expect_equal(unname(z$bands[[sc]][2L, ]), unname(z$bands[[sc]][4L, ]))
  }
})

test_that("the compact response detects a substitution at either fitted scale", {
  old <- options(np.largeh=TRUE, np.largeh.rel.tol=.001)
  on.exit(options(old))
  d <- data.frame(x=seq(-.5,.5,length.out=27), y=seq(-.5,.5,length.out=27))
  calc <- function(h) {
    bw <- npcdistbw(y ~ x, data=d, regtype="lc", cykertype="epanechnikov",
                    bws=unname(rev(h)), bandwidth.compute=FALSE)
    .np_conditional_cdf_centers(bw, d["x"], d["y"], d["x"], d["y"])
  }
  for(axis in 1:2) {
    radius <- if(axis==1) sqrt(-2*log1p(-.001)) else sqrt(.001)
    for(scale in c(1.01, .80)) {
      h <- c(.3,.3)
      h[axis] <- scale/radius
      expect_error(calc(h), "near-constant kernel substitution")
      options(np.largeh=FALSE)
      expect_true(all(is.finite(calc(h))))
      options(np.largeh=TRUE)
    }
  }
})

test_that("compact-response extrema reject signed conditioning weights", {
  z <- list(w=matrix(c(-.1, 1.1)), jw=matrix(c(.2, .8)), ix=1L)
  expect_error(.np_conditional_pairs_range(z, "epanechnikov"), "nonnegative conditioning weights")
})

test_that("the uniform response uses both finite knots of the signed donor", {
  w <- c(0, .03, .6, .37, 0)
  w2 <- c(.07, .5, .2, 0, .23)
  z <- list(w=matrix(w), jw=matrix(w2), ix=1L)
  s <- .np_conditional_pairs_range(z, "uniform")
  knots <- c(-sqrt(2), -1, 1, sqrt(2))
  G <- function(t) pmax(0, pmin(1, (1+t)/2))
  v <- outer(2*w, G(knots)) - outer(w2, G(knots/sqrt(2)))
  lo <- apply(v, 1L, min)
  hi <- apply(v, 1L, max)
  expect_equal(s$lower, sum(lo), tolerance=1e-14)
  expect_equal(drop(s$span*s$q), hi-lo, tolerance=1e-14)
  expect_error(.np_conditional_pairs_range(list(w=matrix(-1),jw=matrix(1),ix=1L),
                                            "uniform"), "nonnegative")
})

test_that("uniform response pairs preserve the fit and reject the widened global limit", {
  d <- data.frame(x=seq(-.5,.5,length.out=51))
  d$y <- .5*cos(seq_along(d$x))
  fit <- function(h) npcdistbw(y~x, data=d, regtype="lc", cykertype="uniform",
                                bws=c(h,.2), bandwidth.compute=FALSE)
  old <- options(np.largeh=TRUE)
  on.exit(options(old))
  args <- list(x=npcdist(fit(.25)), center="bias-corrected", band="all",
               errors="bootstrap", interval="percentile", bootstrap="inid",
               bias.control=np_conditional_bias_control("pairs-score"),
               neval=3L, B=399L, output="data", random.seed=309L)
  set.seed(119L)
  before <- .Random.seed
  v <- do.call(plot,args)
  expect_identical(.Random.seed,before)
  a <- as.data.frame(v,panel=1L,target=1L,band="simultaneous")
  reference <- as.data.frame(plot(args$x,errors="none",neval=3L,output="data"),
                             panel=1L,target=1L)
  expect_equal(a$estimate,reference$estimate,tolerance=2e-8)
  expect_true(all(a$lower <= a$upper))
  bw <- fit(.8*diff(range(d$y)))
  expect_error(.np_conditional_cdf_centers(bw,d["x"],d["y"],d["x"],d["y"]),
               "near-constant")
})

test_that("the local pointwise choice preserves compact response kernels", {
  x <- seq(-1, 1, length.out = 53L)
  d <- data.frame(x = x, y = x + cos(seq_along(x)))
  for (kernel in c("uniform", "epanechnikov")) {
    bw <- npcdistbw(y ~ x, data = d, cykertype = kernel,
                    bws = c(.35, .4), bandwidth.compute = FALSE, regtype = "lc")
    args <- list(x = npcdist(bw), errors = "bootstrap", bootstrap = "inid",
                 interval = "percentile", band = "all", center = "bias-corrected",
                 B = 399L, random.seed = 71L, neval = 3L, output = "data")
    args$bias.control <- np_conditional_bias_control("pairs-score")
    sharp <- do.call(plot, args)
    args$bias.control <- np_conditional_bias_control("pairs-score", pointwise = "local-wilson")
    local <- do.call(plot, args)
    expect_identical(sharp$inference[[1L]]$critical_value_raw,
                     local$inference[[1L]]$critical_value_raw)
    expect_match(local$inference[[1L]]$endpoint_construction, "Wilson")
    point <- as.data.frame(local, panel = 1L, target = 1L, band = "pointwise")
    joint <- as.data.frame(local, panel = 1L, target = 1L, band = "simultaneous")
    expect_true(all(joint$lower <= point$lower & joint$upper >= point$upper))
  }
})

test_that("local Wilson companion uses positive probabilities and preserves the original control", {
  old <- np_conditional_bias_control("pairs-score")
  expect_identical(names(old), c("construction", "version"))
  expect_identical(old$version, 1L)
  new <- np_conditional_bias_control("pairs-score", pointwise = "local-wilson")
  expect_identical(new$version, 2L)
  expect_identical(new, unserialize(serialize(new, NULL)))
  expect_error(np_conditional_bias_control("analytic", pointwise = "local-wilson"), "pairs-score")
  expect_error(np_conditional_bias_control("pairs-score", pointwise = "other"), "local-wilson")
  y <- c(1, 0, 1, 3)
  x <- c(-1, 0, 1, 2)
  grid <- data.frame(x = c(0, 1, 0, 1, 0), y = c(-1, 0, 1, 4, 1))
  got <- .np_conditional_local_wilson(x, y, grid, .8, .05)
  reference <- t(vapply(seq_len(nrow(grid)), function(i) {
    a <- dnorm((x - grid$x[i]) / .8)
    a <- a / sum(a)
    phat <- sum(a[y <= grid$y[i]])
    z2 <- qnorm(.975)^2 * sum(a * a)
    disc <- z2 * phat * (1 - phat) + z2^2 / 4
    (phat + z2 / 2 + c(-1, 1) * sqrt(disc)) / (1 + z2)
  }, numeric(2L)))
  expect_equal(got, reference, tolerance = 1e-14)
  expect_identical(got[3, ], got[5, ])
})
test_that("beta pairs intervals respect response support and the native CDF", {
  x <- seq(-1, 1, length.out = 57L)
  y <- (-5 + 8 * plogis(x + sin(3 * x)))
  d <- data.frame(x = x, y = y)
  g <- data.frame(x = c(-.2, 0, .3, .3), y = c(-5, -1, 3, 3))
  h <- c(.36, .74)
  provider <- .np_conditional_beta_pairs_provider(y, g$y, h[2], c(-5, 3))
  a <- provider(x, y, g, h, .05, "beta")
  expect_equal(as.numeric(a$score[1, ]), c(0, 0))
  expect_equal(as.numeric(a$score[3, ]), c(1, 1))
  expect_identical(a$score[3, ], a$score[4, ])
  expect_error(.np_conditional_beta_pairs_provider(y, -6, h[2], c(-5, 3)), "inside")
  for (scale in c(1, sqrt(2))) {
    bw <- npcdistbw(y ~ x, data = d, regtype = "lc", cykertype = "beta",
                    cykerbound = "fixed", cykerlb = -5, cykerub = 3,
                    bws = rev(h) * scale, bandwidth.compute = FALSE)
    native <- fitted(npcdist(bw, txdat = d["x"], tydat = d["y"],
                            exdat = g["x"], eydat = g["y"]))
    observed <- colSums(if (scale == 1) a$plan$main else a$plan$wide)
    expect_equal(unname(observed), as.numeric(native), tolerance = 2e-12)
  }
  args <- list(x = npcdist(bw), errors = "bootstrap", bootstrap = "inid",
               band = "all", center = "bias-corrected", interval = "percentile",
               B = 399L, neval = 3L, output = "data", random.seed = 17L)
  args$bias.control <- np_conditional_bias_control("pairs-score")
  sharp <- do.call(plot, args)
  args$bias.control <- np_conditional_bias_control("pairs-score", pointwise = "local-wilson")
  local <- do.call(plot, args)
  expect_identical(local$inference[[1L]]$critical_value_raw, sharp$inference[[1L]]$critical_value_raw)
  joint <- as.data.frame(local, panel = 1L, target = 1L, band = "simultaneous")
  point <- as.data.frame(local, panel = 1L, target = 1L, band = "pointwise")
  expect_true(all(joint$lower <= point$lower & joint$upper >= point$upper))
})

test_that("higher-order conditional CDF pairs retain both fitted kernels and dilations", {
  x <- seq(-1.2, 1.5, length.out = 61L)
  d <- data.frame(x = x, y = .5 * x + cos(seq_along(x)))
  grid <- data.frame(x = c(-.3, -.1, .4, -.3), y = c(-.8, .2, .65, -.8))
  h <- c(.36, .41)
  for (kernel in c("gaussian", "epanechnikov")) {
    for (order in c(4L, 6L, 8L)) {
      make <- function(scale) {
        npcdistbw(y ~ x, data = d, cykertype = kernel, cykerorder = order,
                  bws = rev(h * scale), bandwidth.compute = FALSE, regtype = "lc")
      }
      evaluate <- function(b) fitted(npcdist(b, txdat = d["x"], tydat = d["y"],
                                              exdat = grid["x"], eydat = grid["y"]))
      main <- make(1)
      wide <- make(c(sqrt(2), 2^(1 / order)))
      provider <- .np_conditional_high_pairs_provider(d$y, grid$y, h[2L], kernel, order)
      z <- provider(x, d$y, grid, h, .05, kernel)
      expect_equal(colSums(z$plan$main), as.numeric(evaluate(main)), tolerance = 2e-12)
      expect_equal(colSums(z$plan$wide), as.numeric(evaluate(wide)), tolerance = 2e-12)
      expect_identical(z$score[1L, ], z$score[4L, ])
      args <- list(x = main, errors = "bootstrap", bootstrap = "inid",
                   interval = "percentile", center = "bias-corrected", band = "all",
                   B = 99L, neval = 3L, random.seed = 71L, output = "data",
                   bias.control = np_conditional_bias_control("pairs-score", pointwise = "local-wilson"))
      set.seed(937L)
      seed <- .Random.seed
      pd <- suppressWarnings(do.call(plot, args))
      expect_identical(.Random.seed, seed)
      point <- as.data.frame(pd, panel = 1L, target = 1L, band = "pointwise")
      joint <- as.data.frame(pd, panel = 1L, target = 1L, band = "simultaneous")
      expect_true(all(is.finite(point$lower) & joint$lower <= point$lower & point$upper <= joint$upper))
      wrong <- args
      wrong$bootstrap <- "fixed"
      expect_error(do.call(plot, wrong), "inid|unavailable")
    }
  }
})

test_that("bounded beta conditioning retains the fitted rows and all pairs counts", {
  d <- data.frame(x = seq(-3, 2, length.out = 75L), y = seq(.01, .99, length.out = 75L))
  g <- data.frame(x = c(-3, .3, 2, .3), y = c(.05, .4, .95, .4))
  h <- .45
  cache <- .np_conditional_beta_conditioning(d$x, unique(g$x), h, c(-3, 2))
  expect_equal(colSums(cache$main), rep(1, 3L), tolerance = 1e-14)
  q <- (d$x + 3)/5
  want <- vapply((unique(g$x)+3)/5, function(t) {
    w <- dbeta(q, 1 + t*(5/h)^2, 1+(1-t)*(5/h)^2)
    w/sum(w)
  }, numeric(nrow(d)))
  expect_equal(unname(cache$main), unname(want), tolerance = 2e-12)
  expect_error(.np_conditional_beta_conditioning(d$x, 2.1, h, c(-3, 2)), "inside known")
  sparse <- list(targets = 0, main = matrix(c(1, 0)), wide = matrix(c(.5, .5)))
  tasks <- list(list(value = matrix(c(0L, 2L), ncol = 1L)))
  expect_error(.np_conditional_pairs_support(c(0, 1), 0, h, tasks, sparse), "zero or nonfinite")
  for (kernel in c("beta", "gaussian")) {
    args <- list(formula = y ~ x, data = d, regtype = "lc", cxkertype = "beta",
                 cxkerbound = "fixed", cxkerlb = -3, cxkerub = 2,
                 bws = c(.14, h), bandwidth.compute = FALSE)
    if (kernel == "beta")
      args <- c(args, list(cykertype = "beta", cykerbound = "fixed", cykerlb = 0, cykerub = 1))
    bw <- do.call(npcdistbw, args)
    main <- as.numeric(fitted(npcdist(bw, txdat = d["x"], tydat = d["y"],
                                      exdat = g["x"], eydat = g["y"])))
    args$bws <- sqrt(2) * args$bws
    wide <- as.numeric(fitted(npcdist(do.call(npcdistbw, args), txdat = d["x"], tydat = d["y"],
                                      exdat = g["x"], eydat = g["y"])))
    provider <- if (kernel == "beta") .np_conditional_beta_pairs_provider(d$y, g$y, .14, c(0, 1), cache)
                else .np_conditional_high_pairs_provider(d$y, g$y, .14, "gaussian", 2L, cache)
    tile <- provider(d$x, d$y, g, c(h, .14), .05, kernel)
    expect_equal(colSums(tile$plan$main), main, tolerance = 2e-8)
    expect_equal(colSums(tile$plan$wide), wide, tolerance = 2e-8)
    controls <- list(x = bw, errors = "bootstrap", bootstrap = "inid", center = "bias-corrected",
                     band = "all", interval = "percentile", output = "data", neval = 4,
                     B = 99L, random.seed = 78L, bias.control = np_conditional_bias_control("pairs-score", "local-wilson"))
    set.seed(12L)
    before <- .Random.seed
    pd <- suppressWarnings(do.call(plot, controls))
    expect_identical(.Random.seed, before)
    pt <- as.data.frame(pd, panel = 1L, target = 1L, band = "pointwise")
    sm <- as.data.frame(pd, panel = 1L, target = 1L, band = "simultaneous")
    expect_true(all(is.finite(pt$lower) & sm$lower <= pt$lower & pt$upper <= sm$upper))
    controls$bootstrap <- "fixed"
    expect_error(do.call(plot, controls), "inid|unavailable")
  }
})

test_that("conditional pairs use the complete fitted predictor profile", {
  d <- data.frame(x = seq(-1, 1, length.out = 72L),
                  z = ordered(rep(c("-2","0","3"), 24L), levels = c("-2","0","3")))
  d$y <- sin(2 * d$x) + .2 * as.integer(d$z) + .1 * cos(seq_len(nrow(d)))
  b <- npcdistbw(y ~ z + x, data = d, regtype = "lc", oxkertype = "liracine",
                 bws = c(.3, .15, .34), bandwidth.compute = FALSE)
  ex <- d[c(1L, 30L, 48L, 1L), c("z", "x")]
  ex$x[3L] <- ex$x[2L]
  ex$z[3L] <- d$z[1L]
  profiles <- .np_conditional_profile_rows(ex)
  expect_identical(profiles$index, c(1L, 2L, 3L, 1L))
  cache <- .np_conditional_profile_weights(b, d[c("z","x")], ex[profiles$first, , drop = FALSE])
  expect_equal(colSums(cache$main), rep(1, 3L), tolerance = 2e-14, ignore_attr = TRUE)
  expect_equal(colSums(cache$wide), rep(1, 3L), tolerance = 2e-14, ignore_attr = TRUE)
  args <- list(x = npcdist(b), errors = "bootstrap", bootstrap = "inid", interval = "percentile",
               center = "bias-corrected", band = "all", output = "data", neval = 7L, B = 99L,
               random.seed = 55L, bias.control = np_conditional_bias_control("pairs-score","local-wilson"))
  set.seed(29L)
  before <- .Random.seed
  pd <- suppressWarnings(do.call(plot, args))
  expect_identical(.Random.seed, before)
  expect_length(pd$panels, 3L)
  for (j in seq_len(3L)) {
    a <- as.data.frame(pd, panel=j, target=1L, band="pointwise")
    s <- as.data.frame(pd, panel=j, target=1L, band="simultaneous")
    expect_true(all(is.finite(a$lower) & a$lower <= a$upper & s$lower <= a$lower & s$upper >= a$upper))
  }
  wide <- npcdistbw(y ~ z + x, data = d, bws = c(.3, .55, .34), bandwidth.compute = FALSE)
  expect_error(.np_conditional_profile_weights(wide, d[c("z","x")], ex), "outside")
})

test_that("special conditioning weights cannot enter the Gaussian-X-only preparation", {
  y <- seq(-1, 1, length.out = 12L)
  g <- data.frame(x = 1L, y = 0)
  weights <- function() list(targets = 1L,
    main = matrix(rep(1 / length(y), length(y))), wide = matrix(rep(1 / length(y), length(y))))
  set.seed(319L)
  before <- .Random.seed
  expect_error(.np_conditional_pairs_score(y, y, g, c(.3,.2), 99L, .05, NULL, 0,
    response.kernel = "uniform", profile.weights = weights), "no qualified response range")
  expect_identical(.Random.seed, before)
})

test_that("a profile RMS envelope preserves both pointwise choices and the existing band", {
  expect_identical(np_conditional_bias_control("pairs-score")$version, 1L)
  expect_identical(np_conditional_bias_control("pairs-score", "local-wilson")$version, 2L)
  expect_error(np_conditional_bias_control("analytic", simultaneous = "rms-envelope"), "requires construction")
  expect_error(np_conditional_bias_control("pairs-score", simultaneous = "other"), "rms-envelope")
  x <- seq(-1, 1, length.out = 30L)
  dat <- data.frame(x = x, level = factor(rep(letters[1:3], 10L)), y = sin(x) + rep(c(-.4,.3,.1), 10L))
  bw <- npcdistbw(y ~ x + level, data = dat, regtype = "lc",
                  bws = c(.3,.4,.14), bandwidth.compute = FALSE)
  one <- npcdistbw(y ~ x, data = dat, regtype = "lc",
                   bws = c(.3,.4), bandwidth.compute = FALSE)
  args <- list(x = bw, neval = 7L, center = "bias-corrected", errors = "bootstrap",
               interval = "percentile", bootstrap = "inid", band = "all",
               xtrim = .1, B = 199L, random.seed = 51L, output = "data")
  set.seed(27L)
  before <- .Random.seed
  for (point in list(NULL, "local-wilson")) {
    args$bias.control <- np_conditional_bias_control("pairs-score", point)
    plain <- suppressWarnings(do.call(plot, args))
    args$bias.control <- np_conditional_bias_control("pairs-score", point, "rms-envelope")
    augmented <- suppressWarnings(do.call(plot, args))
    for (j in seq_along(plain$panels)) {
      bounds <- function(obj, band) as.data.frame(obj, panel = j, target = 1L, band = band)
      expect_identical(bounds(plain, "pointwise"), bounds(augmented, "pointwise"))
      a <- bounds(plain, "simultaneous")
      b <- bounds(augmented, "simultaneous")
      c <- bounds(augmented, "bonferroni")
      expect_true(all(c$lower <= b$lower & b$lower <= a$lower & a$upper <= b$upper & b$upper <= c$upper))
      components <- augmented$inference[[j]]$calibration_details$critical_components
      expect_length(components$root.rms, nrow(b))
      expect_true(is.finite(components$rms.multiplier) && components$rms.multiplier >= 0)
    }
  }
  args$x <- one
  expect_error(do.call(plot, args), "two-predictor")
  expect_identical(.Random.seed, before)
})

test_that("an RMS request exceeding its memory allowance refuses before preparing weights", {
  grid <- data.frame(x = c(0, 1), y = c(0, 1))
  fail <- function() stop("unexpected call to prepare profile weights")
  set.seed(83L)
  before <- .Random.seed
  expect_error(.np_conditional_pairs_score(1:8, 1:8, grid, c(.3,.4), 99L, .05,
    NULL, 256 * 1024^2 - 1, profile.weights = fail, simultaneous = "rms-envelope"), "RMS envelope exceeds")
  expect_identical(.Random.seed, before)
})

test_that("RMS pairs preserve changes below one CDF rounding unit", {
  response <- c(0, 0, 0, .5)
  target <- seq(10.5, 12.5, by=.1)
  h <- c(.3, 1)
  w <- matrix(c(.125, .375, .25, .25), ncol=1L)
  v <- matrix(c(.25, .25, .125, .375), ncol=1L)
  weights <- list(targets=0, main=w, wide=v)
  prepare <- .np_conditional_high_pairs_provider(response, target, 1, "gaussian",
                                                  2L, weights, rms=TRUE)
  plan <- prepare(1:4, response, data.frame(x=0, y=target), h, .05, "gaussian")$plan
  counts <- cbind(rep(1L, 4L), c(2L,1L,1L,0L), c(0L,0L,1L,3L))
  both <- .np_conditional_pairs_block(list(value=counts),plan)
  main <- .np_conditional_high_integral(response, target, 1, "gaussian", 2L)
  wide <- .np_conditional_high_integral(response, target, sqrt(2), "gaussian", 2L)
  change <- function(weights) {
    counts[4L,] * weights[4L] / as.numeric(crossprod(counts,weights)) - weights[4L]
  }
  # In a two-value donor population the exact change factors into a change
  # of normalized mass and a difference of the two stored probabilities.
  want <- 2 * outer(change(w),main[4L,]-main[1L,]) -
              outer(change(v),wide[4L,]-wide[1L,])
  m <- length(target)
  expect_equal(dim(both), c(3L, 2L*m))
  stable <- both[,m+seq_len(m),drop=FALSE]
  tiny <- abs(want)>0 & abs(want)<.Machine$double.eps
  expect_true(any(tiny))
  expect_true(all(abs(stable[tiny]/want[tiny]-1)<1e-12))
  expect_true(all(abs(stable-want)<=1e-12*abs(want)))
  plan$stable <- NULL
  expect_identical(both[,seq_len(m),drop=FALSE],
                   .np_conditional_pairs_block(list(value=counts),plan))
})
