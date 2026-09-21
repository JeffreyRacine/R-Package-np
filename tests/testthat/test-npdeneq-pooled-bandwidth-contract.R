dh_payload <- function(x) x[c("In", "Tn", "In.bootstrap", "Tn.bootstrap", "In.P", "Tn.P")]
dh_quiet <- function() {
  old <- options(np.messages = FALSE)
  withr::defer(options(old), envir = parent.frame())
}

test_that("reference adjustment retains compact symmetric kernel choices", {

  dh_quiet()
  d <- data.frame(v = seq_len(30)/30)
  adjust <- getFromNamespace(".npdeneq_reference_bandwidth", "np")
  for (kernel in c("gaussian", "epanechnikov", "uniform")) {
    bw <- npudensbw(dat = d, bws = .2, bandwidth.compute = FALSE, ckertype = kernel)
    if (kernel == "uniform")
      expect_warning(actual <- adjust(bw, 10, 20), "ignoring kernel order")
    else
      expect_warning(actual <- adjust(bw, 10, 20), NA)
    expect_identical(actual$ckertype, kernel)
    expect_equal(unname(actual$bw), .2*(30/(2*10*20/30))^.2, tolerance = 1e-14)
  }
})

test_that("harmonic reference scaling retains pooled spread and kernel metadata", {
  dh_quiet()
  d <- data.frame(v = seq(-1, 1, length.out = 30),
    w = sin(seq_len(30)), u = factor(rep(1:3, 10)),
    o = ordered(rep(1:3, 10), levels = 1:4))
  adjust <- getFromNamespace(".npdeneq_reference_bandwidth", "np")
  kbw <- getFromNamespace("kbandwidth", "np")
  for (order in c(2, 4, 6, 8)) for (cols in list("v", c("v", "w"),
      c("u", "o"), c("v", "u", "o"))) {
    z <- d[cols]
    raw <- ifelse(vapply(z, is.factor, logical(1)), .1, .3)
    for (scaled in c(FALSE, TRUE)) {
      bw <- npudensbw(dat = z, bws = raw, bwscaling = scaled,
        bandwidth.compute = FALSE, ckerorder = order)
      before <- serialize(bw, NULL)
      base <- kbw(bw)
      actual <- adjust(bw, 10L, 20L)
      # Independent reference-size conversion of the retained physical values.
      H <- 2 * 10 * 20 / 30
      q <- sum(!vapply(z, is.factor, logical(1)))
      powers <- ifelse(vapply(z, is.factor, logical(1)), 2, 1)/(2*order+q)
      expected <- base$bw * (30/H)^powers
      expect_equal(actual$bw, expected, tolerance = 1e-14)
      expect_identical(actual$bw, adjust(bw, 20L, 10L)$bw)
      for (name in c("xmcv", "xdati", "icon", "iord", "iuno",
                     "ckertype", "ckerorder", "ukertype", "okertype", "nobs"))
        expect_identical(actual[[name]], base[[name]])
      expect_false(actual$scaling)
      expect_identical(serialize(bw, NULL), before)
    }
  }
})

test_that("categorical reference scaling explicitly respects kernel constraints", {
  dh_quiet()
  d <- data.frame(u = factor(rep(1:3, 10)), o = ordered(rep(1:3, 10)))
  adjust <- getFromNamespace(".npdeneq_reference_bandwidth", "np")
  for (u in c("aitchisonaitken", "liracine"))
    for (o in c("liracine", "wangvanryzin")) {
      upper <- c(if (u == "aitchisonaitken") 2/3 else 1, 1)
      bw <- npudensbw(dat = d, bws = upper * .9, bandwidth.compute = FALSE,
                     ukertype = u, okertype = o)
      before <- serialize(bw, NULL)
      expect_warning(actual <- adjust(bw, 10, 20),
                     "categorical kernel upper bound for: u, o", fixed = TRUE)
      expect_equal(unname(actual$bw), upper, tolerance = 0)
      expect_identical(serialize(bw, NULL), before)
      zero <- npudensbw(dat = d, bws = c(0, 0), bandwidth.compute = FALSE,
                       ukertype = u, okertype = o)
      expect_warning(out <- adjust(zero, 10, 20), NA)
      expect_identical(unname(out$bw), c(0, 0))
      # Below/at/above the adjustment boundary; no blanket clipping.
      factor <- sqrt(30/(2*10*20/30))
      for (mult in c(1 - 1e-8, 1, 1 + 1e-8)) {
        edge <- npudensbw(dat = d, bws = upper/factor * mult,
                         bandwidth.compute = FALSE, ukertype = u, okertype = o)
        if (mult > 1)
          expect_warning(out <- adjust(edge, 10, 20), "categorical kernel upper bound")
        else
          expect_warning(out <- adjust(edge, 10, 20), NA)
        expect_equal(unname(out$bw), pmin(upper * mult, upper), tolerance = 1e-14)
      }
    }
})

test_that("automatic selection sees both complete samples exactly once", {
  dh_quiet()
  x <- data.frame(v = seq_len(12)/12, u = factor(rep(c("a", "b"), 6)))
  y <- data.frame(v = seq_len(20)/20, u = factor(rep(c("a", "c"), 10)))
  x$v[1] <- NA; y$v[2] <- NA
  pool <- rbind(stats::na.omit(x), stats::na.omit(y))
  state <- new.env(parent = emptyenv()); state$calls <- list()
  selected <- npudensbw(dat = pool, bws = c(.3, .1), bandwidth.compute = FALSE)
  local_mocked_bindings(.npdeneq_select_bandwidth = function(x, ...) {
    state$calls <- c(state$calls, list(x)); selected
  }, .package = "np")
  set.seed(42); before <- .Random.seed
  actual <- npdeneqtest(x, y, B = 9)
  expect_length(state$calls, 1L)
  expect_equal(state$calls[[1L]]$v, pool$v)
  expect_identical(levels(state$calls[[1L]]$u), c("a", "b", "c"))
  expect_identical(.Random.seed, before)
  ratio <- 30/(2*11*19/30)
  control <- npudensbw(dat = pool, bws = c(.3*ratio^.2, .1*ratio^.4),
                      bandwidth.compute = FALSE)
  expected <- npdeneqtest(stats::na.omit(x), stats::na.omit(y), bw.x = control, B = 9)
  expect_equal(dh_payload(actual), dh_payload(expected), tolerance = 2e-12)
})

test_that("actual pooled search agrees with an explicitly adjusted common bandwidth", {
  dh_quiet()
  x <- data.frame(v = sin(seq_len(18)))
  y <- data.frame(v = cos(seq_len(30)) + .2)
  pool <- rbind(x, y)
  for (scaled in c(FALSE, TRUE)) {
    set.seed(719)
    pooled <- npudensbw(dat = pool, bwmethod = "cv.ls", nmulti = 1,
                       bwscaling = scaled)
    physical <- unlist(pooled$bandwidth, use.names = FALSE)
    reference <- npudensbw(dat = pool, bandwidth.compute = FALSE,
      bws = physical * (48/(2*18*30/48))^.2)
    expected <- npdeneqtest(x, y, bw.x = reference, B = 9)
    expected.seed <- .Random.seed
    set.seed(719)
    actual <- npdeneqtest(x, y, B = 9, nmulti = 1, bwscaling = scaled)
    expect_equal(dh_payload(actual), dh_payload(expected), tolerance = 2e-12)
    expect_identical(.Random.seed, expected.seed)
    expect_true(all(is.finite(c(actual$In, actual$Tn.bootstrap))))
  }
})

test_that("manual scaled bandwidths retain their original first-sample units", {
  dh_quiet()
  x <- data.frame(v = seq_len(12)/12)
  y <- data.frame(v = seq_len(28)/7)
  bw <- npudensbw(dat = x, bws = .8, bandwidth.compute = FALSE, bwscaling = TRUE)
  expected <- npdeneqtest(x, y, bw.x = bw, B = 9)
  actual <- npdeneqtest(x, y, bws = .8, bandwidth.compute = FALSE,
                       bwscaling = TRUE, B = 9)
  expect_identical(dh_payload(actual), dh_payload(expected))
  for (manual in list(0, "FALSE")) {
    actual <- npdeneqtest(x, y, bws = .8, bandwidth.compute = manual,
                         bwscaling = TRUE, B = 9)
    expect_identical(dh_payload(actual), dh_payload(expected))
  }
  actual <- npdeneqtest(x, y, bws = .8, bandwidth.comp = FALSE,
                       bwscaling = TRUE, B = 9)
  expect_identical(dh_payload(actual), dh_payload(expected))
  local_mocked_bindings(.npdeneq_select_bandwidth = function(...) stop("SEARCH_ENTERED"),
                        .package = "np")
  expect_identical(dh_payload(npdeneqtest(x, y, bw.y = bw, B = 9)),
                   dh_payload(expected))
})
