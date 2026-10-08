library(np)

collect_plot_payload_fields <- function(x, path = character()) {
  out <- list()

  if (is.null(x)) {
    return(out)
  }

  leaf <- if (length(path)) path[[length(path)]] else ""
  if (is.data.frame(x) && leaf %in% c("eval", "xeval", "yeval", "evalx", "evaly", "evalz")) {
    out[[paste(path, collapse = ".")]] <- x
    return(out)
  }
  if (is.data.frame(x)) x <- as.list(x)

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
  if (!(leaf %in% c("mean", "merr", "dens", "dist", "condens", "condist",
                    "derr", "conderr", "eval", "xeval", "yeval",
                    "evalx", "evaly", "evalz"))) {
    return(out)
  }

  out[[paste(path, collapse = ".")]] <-
    if (leaf %in% c("derr", "conderr")) x else as.numeric(x)
  out
}

run_exact_frozen_plot_pair <- function(fit, ..., seed = 9001L) {
  extra <- list(...)
  exact <- NULL
  frozen <- NULL

  set.seed(seed)
  suppressWarnings(capture.output(
    exact <- do.call(plot, c(list(fit), extra, list(boot.control = np_boot_control(nonfixed = "exact"))))
  ))

  set.seed(seed)
  suppressWarnings(capture.output(
    frozen <- do.call(plot, c(list(fit), extra, list(boot.control = np_boot_control(nonfixed = "frozen"))))
  ))

  list(exact = exact, frozen = frozen)
}

expect_plot_payload_comparable <- function(pair,
                                           label,
                                           min.merr.corr,
                                           max.merr.rel) {
  exact.fields <- collect_plot_payload_fields(pair$exact)
  frozen.fields <- collect_plot_payload_fields(pair$frozen)

  expect_true(length(exact.fields) > 0L, info = label)
  expect_true(length(frozen.fields) > 0L, info = label)
  expect_equal(sort(names(exact.fields)), sort(names(frozen.fields)), info = label)

  required <- list(npudens = c("dens", "eval", "derr"),
                   npudist = c("dist", "eval", "derr"),
                   npcdens = c("condens", "xeval", "yeval", "conderr"),
                   npcdist = c("condist", "xeval", "yeval", "conderr"))
  for (field in required[[label]])
    expect_true(any(grepl(paste0("(^|[.])", field, "$"), names(exact.fields))),
                info = paste(label, "missing", field))

  if (identical(label, "npplreg")) {
    expect_true(any(grepl("(^|[.])evalx$", names(exact.fields))), info = label)
    expect_true(any(grepl("(^|[.])evalz$", names(exact.fields))), info = label)
  }
  for (nm in names(exact.fields)) {
    exact.val <- exact.fields[[nm]]
    frozen.val <- frozen.fields[[nm]]
    leaf <- sub("^.*\\.", "", nm)
    if (is.data.frame(exact.val)) {
      expect_identical(dim(frozen.val), dim(exact.val), info = label)
      expect_identical(names(frozen.val), names(exact.val), info = label)
      expect_true(all(is.finite(as.matrix(exact.val))), info = label)
      expect_true(all(is.finite(as.matrix(frozen.val))), info = label)
      expect_equal(frozen.val, exact.val, tolerance = 1e-10, info = label)
      next
    }

    expect_equal(length(frozen.val), length(exact.val), info = sprintf("%s %s length", label, nm))
    if (leaf %in% c("mean", "dens", "dist", "condens", "condist",
                    "eval", "xeval", "yeval", "evalx", "evaly", "evalz")) {
      expect_true(all(is.finite(exact.val)), info = sprintf("%s %s exact finite", label, nm))
      expect_true(all(is.finite(frozen.val)), info = sprintf("%s %s frozen finite", label, nm))
      expect_equal(frozen.val, exact.val, tolerance = 1e-10, info = sprintf("%s %s", label, nm))
      next
    }

    if (leaf %in% c("derr", "conderr")) {
      # These fields were not collected by the old test. Exact resampling
      # recomputes NN bandwidths; frozen resampling holds them fixed. There
      # is no contract that their uncertainty estimates meet the regression
      # merr correlation thresholds below. Check the interval payload, while
      # the point estimates and coordinates above must agree exactly.
      expect_identical(dim(frozen.val), dim(exact.val), info = label)
      expect_identical(is.na(frozen.val), is.na(exact.val), info = label)
      for (value in list(exact.val, frozen.val)) {
        observed <- value[!is.na(value)]
        expect_true(length(observed) > 0L, info = label)
        expect_true(all(is.finite(observed)), info = label)
        expect_true(any(observed != 0), info = label)
      }
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
    output = "data",
    errors = "bootstrap",
    bootstrap = "inid",
    B = 41L,
    band = "pointwise",
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
    output = "data",
    errors = "bootstrap",
    bootstrap = "inid",
    B = 41L,
    band = "pointwise",
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
    output = "data",
    errors = "bootstrap",
    bootstrap = "geom",
    B = 41L,
    band = "pointwise",
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
    output = "data",
    errors = "bootstrap",
    bootstrap = "inid",
    B = 41L,
    band = "pointwise",
    neval = 8L
  )
  expect_plot_payload_comparable(sc.pair, "npscoef", min.merr.corr = 0.95, max.merr.rel = 0.35)
})

test_that("density and distribution plots preserve points, grids and interval layout", {
  set.seed(20260323)

  n.u <- 75L
  y.u <- runif(n.u, -1, 1) + rnorm(n.u, sd = 0.1)
  ud.fit <- npudens(~ y.u, nmulti = 1, bwtype = "adaptive_nn")
  ud.pair <- run_exact_frozen_plot_pair(
    ud.fit,
    output = "data",
    errors = "bootstrap",
    bootstrap = "inid",
    B = 41L,
    band = "pointwise",
    neval = 8L
  )
  expect_plot_payload_comparable(ud.pair, "npudens", min.merr.corr = 0.95, max.merr.rel = 0.75)

  uf.fit <- npudist(~ y.u, nmulti = 1, bwtype = "adaptive_nn")
  uf.pair <- run_exact_frozen_plot_pair(
    uf.fit,
    output = "data",
    errors = "bootstrap",
    bootstrap = "inid",
    B = 41L,
    band = "pointwise",
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
    output = "data",
    errors = "bootstrap",
    bootstrap = "inid",
    B = 41L,
    band = "pointwise",
    neval = 8L
  )
  expect_plot_payload_comparable(cd.pair, "npcdens", min.merr.corr = 0.55, max.merr.rel = 1.20)

  cf.fit <- npcdist(y ~ x, data = c.dat, nmulti = 1, bwtype = "generalized_nn")
  cf.pair <- run_exact_frozen_plot_pair(
    cf.fit,
    view = "fixed",
    output = "data",
    errors = "bootstrap",
    bootstrap = "inid",
    B = 41L,
    band = "pointwise",
    neval = 8L
  )
  expect_plot_payload_comparable(cf.pair, "npcdist", min.merr.corr = 0.90, max.merr.rel = 0.75)
})

test_that("density plot interval values equal quantiles of their own bootstrap draws", {
  pkg<-getNamespaceName(environment(npudens))
  original<-getFromNamespace('.np_plot_bootstrap_centered_interval_payload',pkg)
  seen<-list()
  testthat::local_mocked_bindings(.np_plot_bootstrap_centered_interval_payload=
    function(boot.t,t0,alpha,band.type,center,...) {
      expect_identical(band.type,'pointwise');expect_identical(center,'estimate')
      expect_identical(alpha, .05)
      expect_equal(nrow(as.matrix(boot.t)), 41L)
      bounds<-t(apply(as.matrix(boot.t),2,quantile,probs=c(alpha/2,1-alpha/2)))
      seen[[length(seen)+1L]]<<-list(point=as.numeric(t0),err=sweep(bounds,1,as.numeric(t0),'-'))
      original(boot.t=boot.t,t0=t0,alpha=alpha,band.type=band.type,center=center,...)
    },.package=pkg)
  set.seed(36174);n<-80L;d<-data.frame(x=runif(n));d$y<-sin(3*d$x)+rnorm(n,sd=.4)
  fits<-list(npudens(npudensbw(~y,data=d,bws=20,bwtype='adaptive_nn',bandwidth.compute=FALSE)),
    npudist(npudistbw(~y,data=d,bws=20,bwtype='adaptive_nn',bandwidth.compute=FALSE)),
    npcdens(npcdensbw(y~x,data=d,bws=c(20,20),bwtype='generalized_nn',bandwidth.compute=FALSE)),
    npcdist(npcdistbw(y~x,data=d,bws=c(20,20),bwtype='generalized_nn',bandwidth.compute=FALSE)))
  for (fit in fits) for(mode in c('exact','frozen')) {
    seen<-list();set.seed(9174)
    suppressWarnings(capture.output(out<-plot(fit,output='data',perspective=FALSE,view='fixed',
      errors='bootstrap',bootstrap='inid',B=41L,band='pointwise',center='estimate',neval=7L,
      boot.control=np_boot_control(nonfixed=mode))))
    fields<-collect_plot_payload_fields(out)
    intervals<-grep('(^|[.])(derr|conderr)$',names(fields),value=TRUE)
    expect_true(length(intervals)>0L);expect_equal(length(intervals),length(seen))
    for(nm in intervals){
      prefix<-sub('(derr|conderr)$','',nm)
      point.names<-paste0(prefix,c('dens','dist','condens','condist'))
      point<-fields[[intersect(point.names,names(fields))]]
      index<-which(vapply(seen,function(s)length(s$point)==length(point)&&max(abs(s$point-point))<1e-10,logical(1)))
      expect_length(index,1L)
      expected<-seen[[index]]$err;actual<-fields[[nm]]
      expect_equal(unname(actual),unname(expected),tolerance=1e-12)
      # Calibrate this fixture against the reported undetected mutations.
      expect_false(isTRUE(all.equal(unname(actual*3),unname(expected),tolerance=1e-12)))
      expect_false(isTRUE(all.equal(unname(actual[,2:1,drop=FALSE]),unname(expected),tolerance=1e-12)))
    }
  }
})
