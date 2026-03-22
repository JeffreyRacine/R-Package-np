library(np)

collect_plot_payload_fields <- function(x, path = character()) {
  out <- list()

  if (is.null(x)) {
    return(out)
  }

  if (is.data.frame(x)) {
    x <- as.list(x)
  }

  if (is.list(x)) {
    nm <- names(x)
    if (is.null(nm)) {
      nm <- as.character(seq_along(x))
    }

    for (i in seq_along(x)) {
      out <- c(out, collect_plot_payload_fields(x[[i]], c(path, nm[i])))
    }
    return(out)
  }

  if (!is.numeric(x)) {
    return(out)
  }

  leaf <- if (length(path)) path[[length(path)]] else ""
  if (!(leaf %in% c("mean", "merr", "evalx", "evaly", "evalz"))) {
    return(out)
  }

  out[[paste(path, collapse = ".")]] <- as.numeric(x)
  out
}

run_exact_frozen_plot_pair <- function(fit, ..., seed = 9001L) {
  extra <- list(...)
  exact <- NULL
  frozen <- NULL

  set.seed(seed)
  suppressWarnings(capture.output(
    exact <- do.call(plot, c(list(fit), extra, list(plot.errors.boot.nonfixed = "exact")))
  ))

  set.seed(seed)
  suppressWarnings(capture.output(
    frozen <- do.call(plot, c(list(fit), extra, list(plot.errors.boot.nonfixed = "frozen")))
  ))

  list(exact = exact, frozen = frozen)
}

expect_plot_payload_comparable <- function(pair,
                                           label,
                                           min.merr.corr,
                                           max.merr.rel) {
  exact.fields <- collect_plot_payload_fields(pair$exact)
  frozen.fields <- collect_plot_payload_fields(pair$frozen)

  expect_equal(sort(names(exact.fields)), sort(names(frozen.fields)), info = label)

  for (nm in names(exact.fields)) {
    exact.val <- exact.fields[[nm]]
    frozen.val <- frozen.fields[[nm]]
    leaf <- sub("^.*\\.", "", nm)

    expect_equal(length(frozen.val), length(exact.val), info = sprintf("%s %s length", label, nm))
    if (leaf %in% c("mean", "evalx", "evaly", "evalz")) {
      expect_true(all(is.finite(exact.val)), info = sprintf("%s %s exact finite", label, nm))
      expect_true(all(is.finite(frozen.val)), info = sprintf("%s %s frozen finite", label, nm))
      expect_equal(frozen.val, exact.val, tolerance = 1e-10, info = sprintf("%s %s", label, nm))
      next
    }

    expect_equal(is.na(frozen.val), is.na(exact.val), info = sprintf("%s %s NA mask", label, nm))
    keep <- is.finite(exact.val) & is.finite(frozen.val)
    if (!any(keep)) {
      next
    }

    exact.keep <- exact.val[keep]
    frozen.keep <- frozen.val[keep]
    scale <- max(1e-8, max(abs(exact.keep)))
    rel.max <- max(abs(frozen.keep - exact.keep)) / scale
    corr <- if (length(exact.keep) > 1L &&
                stats::sd(exact.keep) > 0 &&
                stats::sd(frozen.keep) > 0) {
      stats::cor(exact.keep, frozen.keep)
    } else {
      1
    }

    expect_true(corr >= min.merr.corr, info = sprintf("%s %s corr", label, nm))
    expect_true(rel.max <= max.merr.rel, info = sprintf("%s %s rel", label, nm))
  }
}

test_that("exact and frozen plot payloads stay comparable for regression and semiparametric families", {
  set.seed(20260322)

  n.reg <- 70L
  x.reg <- runif(n.reg, -1, 1)
  y.reg <- x.reg + rnorm(n.reg, sd = 0.15)
  reg.fit <- npreg(y.reg ~ x.reg,
                   nmulti = 1,
                   regtype = "lp",
                   degree = 1,
                   bwtype = "adaptive_nn")
  reg.pair <- run_exact_frozen_plot_pair(
    reg.fit,
    view = "fixed",
    gradients = TRUE,
    plot.behavior = "data",
    plot.errors.method = "bootstrap",
    plot.errors.boot.method = "inid",
    plot.errors.boot.num = 41L,
    plot.errors.type = "pointwise",
    neval = 8L
  )
  expect_plot_payload_comparable(reg.pair, "npreg", min.merr.corr = 0.90, max.merr.rel = 0.75)

  n.idx <- 70L
  x.idx <- runif(n.idx, -1, 1)
  z.idx <- rnorm(n.idx)
  y.idx <- x.idx + 0.5 * z.idx + rnorm(n.idx, sd = 0.15)
  idx.fit <- npindex(y.idx ~ x.idx + z.idx,
                     nmulti = 1,
                     gradients = TRUE,
                     bwtype = "adaptive_nn")
  idx.pair <- run_exact_frozen_plot_pair(
    idx.fit,
    view = "fixed",
    plot.behavior = "data",
    plot.errors.method = "bootstrap",
    plot.errors.boot.method = "inid",
    plot.errors.boot.num = 41L,
    plot.errors.type = "pointwise",
    neval = 8L
  )
  expect_plot_payload_comparable(idx.pair, "npindex", min.merr.corr = 0.90, max.merr.rel = 0.80)

  n.pl <- 75L
  x.pl <- runif(n.pl, -1, 1)
  z.pl <- rnorm(n.pl)
  y.pl <- x.pl^2 + z.pl + rnorm(n.pl, sd = 0.2)
  pl.fit <- npplreg(y.pl ~ x.pl | z.pl, nmulti = 1, bwtype = "adaptive_nn")
  pl.pair <- run_exact_frozen_plot_pair(
    pl.fit,
    view = "fixed",
    plot.behavior = "data",
    plot.errors.method = "bootstrap",
    plot.errors.boot.method = "geom",
    plot.errors.boot.num = 41L,
    plot.errors.type = "pointwise",
    neval = 8L
  )
  expect_plot_payload_comparable(pl.pair, "npplreg", min.merr.corr = 0.70, max.merr.rel = 0.55)

  n.sc <- 75L
  x.sc <- runif(n.sc, -1, 1)
  z.sc <- rnorm(n.sc)
  y.sc <- x.sc^2 + z.sc + rnorm(n.sc, sd = 0.2 * stats::sd(x.sc))
  sc.fit <- npscoef(y.sc ~ x.sc | z.sc,
                    nmulti = 1,
                    regtype = "ll",
                    bwtype = "adaptive_nn")
  sc.pair <- run_exact_frozen_plot_pair(
    sc.fit,
    view = "fixed",
    coef = FALSE,
    plot.behavior = "data",
    plot.errors.method = "bootstrap",
    plot.errors.boot.method = "inid",
    plot.errors.boot.num = 41L,
    plot.errors.type = "pointwise",
    neval = 8L
  )
  expect_plot_payload_comparable(sc.pair, "npscoef", min.merr.corr = 0.95, max.merr.rel = 0.35)
})

test_that("exact and frozen plot payloads stay comparable for density and distribution families", {
  set.seed(20260323)

  n.u <- 75L
  y.u <- runif(n.u, -1, 1) + rnorm(n.u, sd = 0.1)
  ud.fit <- npudens(~ y.u, nmulti = 1, bwtype = "adaptive_nn")
  ud.pair <- run_exact_frozen_plot_pair(
    ud.fit,
    plot.behavior = "data",
    plot.errors.method = "bootstrap",
    plot.errors.boot.method = "inid",
    plot.errors.boot.num = 41L,
    plot.errors.type = "pointwise",
    neval = 8L
  )
  expect_plot_payload_comparable(ud.pair, "npudens", min.merr.corr = 0.95, max.merr.rel = 0.75)

  uf.fit <- npudist(~ y.u, nmulti = 1, bwtype = "adaptive_nn")
  uf.pair <- run_exact_frozen_plot_pair(
    uf.fit,
    plot.behavior = "data",
    plot.errors.method = "bootstrap",
    plot.errors.boot.method = "inid",
    plot.errors.boot.num = 41L,
    plot.errors.type = "pointwise",
    neval = 8L
  )
  expect_plot_payload_comparable(uf.pair, "npudist", min.merr.corr = 0.95, max.merr.rel = 0.45)

  n.c <- 75L
  x.c <- runif(n.c, -1, 1)
  y.c <- x.c + rnorm(n.c, sd = 0.2)
  c.dat <- data.frame(x = x.c, y = y.c)
  cd.fit <- npcdens(y ~ x, data = c.dat, nmulti = 1, bwtype = "generalized_nn")
  cd.pair <- run_exact_frozen_plot_pair(
    cd.fit,
    view = "fixed",
    plot.behavior = "data",
    plot.errors.method = "bootstrap",
    plot.errors.boot.method = "inid",
    plot.errors.boot.num = 41L,
    plot.errors.type = "pointwise",
    neval = 8L
  )
  expect_plot_payload_comparable(cd.pair, "npcdens", min.merr.corr = 0.55, max.merr.rel = 1.20)

  cf.fit <- npcdist(y ~ x, data = c.dat, nmulti = 1, bwtype = "generalized_nn")
  cf.pair <- run_exact_frozen_plot_pair(
    cf.fit,
    view = "fixed",
    plot.behavior = "data",
    plot.errors.method = "bootstrap",
    plot.errors.boot.method = "inid",
    plot.errors.boot.num = 41L,
    plot.errors.type = "pointwise",
    neval = 8L
  )
  expect_plot_payload_comparable(cf.pair, "npcdist", min.merr.corr = 0.90, max.merr.rel = 0.75)
})
