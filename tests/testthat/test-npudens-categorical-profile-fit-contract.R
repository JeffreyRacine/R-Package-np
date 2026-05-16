test_that("npudens all-categorical tree profile fit matches dense fit", {
  old.tree <- getOption("np.tree")
  old.compress <- getOption("np.categorical.compress")
  old.messages <- getOption("np.messages")
  on.exit({
    options(np.tree = old.tree)
    options(np.categorical.compress = old.compress)
    options(np.messages = old.messages)
  }, add = TRUE)
  options(np.messages = FALSE)

  run_case <- function(kind, ukertype, okertype, bws, seed) {
    set.seed(seed)
    n <- 960L
    ne <- 240L

    if (identical(kind, "unordered")) {
      dat <- data.frame(
        u1 = factor(sample(letters[1:4], n, TRUE)),
        u2 = factor(sample(LETTERS[1:3], n, TRUE))
      )
      edat <- dat[sample.int(n, ne, TRUE), , drop = FALSE]
      form <- ~ u1 + u2
    } else if (identical(kind, "ordered")) {
      dat <- data.frame(
        o1 = ordered(sample(1:6, n, TRUE)),
        o2 = ordered(sample(1:5, n, TRUE))
      )
      edat <- dat[sample.int(n, ne, TRUE), , drop = FALSE]
      form <- ~ o1 + o2
    } else {
      dat <- data.frame(
        u1 = factor(sample(letters[1:4], n, TRUE)),
        o1 = ordered(sample(1:6, n, TRUE))
      )
      edat <- dat[sample.int(n, ne, TRUE), , drop = FALSE]
      form <- ~ u1 + o1
    }

    options(np.tree = FALSE, np.categorical.compress = FALSE)
    bw.dense <- npudensbw(
      form,
      data = dat,
      edat = edat,
      bws = bws,
      bandwidth.compute = FALSE,
      ukertype = ukertype,
      okertype = okertype
    )
    fit.dense <- npudens(bws = bw.dense)

    options(np.tree = FALSE, np.categorical.compress = TRUE)
    bw.profile <- npudensbw(
      form,
      data = dat,
      edat = edat,
      bws = bws,
      bandwidth.compute = FALSE,
      ukertype = ukertype,
      okertype = okertype
    )
    fit.profile <- npudens(bws = bw.profile)

    expect_equal(
      fitted(fit.profile),
      fitted(fit.dense),
      tolerance = 1e-12,
      info = paste(kind, ukertype, okertype, "density")
    )
    expect_equal(
      se(fit.profile),
      se(fit.dense),
      tolerance = 1e-12,
      info = paste(kind, ukertype, okertype, "se")
    )
  }

  run_case("unordered", "aitchisonaitken", "wangvanryzin", c(0.25, 0.35), 1001L)
  run_case("unordered", "liracine", "wangvanryzin", c(0.25, 0.35), 1002L)
  run_case("ordered", "aitchisonaitken", "wangvanryzin", c(0.25, 0.35), 1003L)
  run_case("ordered", "aitchisonaitken", "liracine", c(0.25, 0.35), 1004L)
  run_case("ordered", "aitchisonaitken", "racineliyan", c(0.25, 0.35), 1005L)
  run_case("mixed", "liracine", "racineliyan", c(0.25, 0.35), 1006L)
})

test_that("one-coordinate categorical tree profiles match dense density and distribution fits", {
  old.tree <- getOption("np.tree")
  old.compress <- getOption("np.categorical.compress")
  old.messages <- getOption("np.messages")
  on.exit({
    options(np.tree = old.tree)
    options(np.categorical.compress = old.compress)
    options(np.messages = old.messages)
  }, add = TRUE)
  options(np.messages = FALSE)

  run_case <- function(kind, bws, seed) {
    set.seed(seed)
    n <- 1000L
    ne <- 120L

    if (identical(kind, "unordered")) {
      dat <- data.frame(x = factor(sample(letters[1:5], n, TRUE)))
      edat <- data.frame(x = factor(letters[1:5], levels = levels(dat$x)))
    } else {
      dat <- data.frame(x = ordered(sample(0:10, n, TRUE)))
      edat <- data.frame(x = ordered(0:10, levels = levels(dat$x)))
    }

    edat <- edat[sample.int(nrow(edat), ne, TRUE), , drop = FALSE]

    options(np.tree = FALSE, np.categorical.compress = FALSE)
    dens.bw.dense <- npudensbw(~ x, data = dat, edat = edat,
                               bws = bws, bandwidth.compute = FALSE)
    dens.dense <- npudens(bws = dens.bw.dense)
    if (identical(kind, "ordered")) {
      dist.bw.dense <- npudistbw(~ x, data = dat, edat = edat,
                                 bws = bws, bandwidth.compute = FALSE)
      dist.dense <- npudist(bws = dist.bw.dense)
    }

    options(np.tree = FALSE, np.categorical.compress = TRUE)
    dens.bw.profile <- npudensbw(~ x, data = dat, edat = edat,
                                 bws = bws, bandwidth.compute = FALSE)
    dens.profile <- npudens(bws = dens.bw.profile)
    if (identical(kind, "ordered")) {
      dist.bw.profile <- npudistbw(~ x, data = dat, edat = edat,
                                   bws = bws, bandwidth.compute = FALSE)
      dist.profile <- npudist(bws = dist.bw.profile)
    }

    expect_equal(fitted(dens.profile), fitted(dens.dense),
                 tolerance = 1e-12, info = paste(kind, "density"))
    expect_equal(se(dens.profile), se(dens.dense),
                 tolerance = 1e-12, info = paste(kind, "density se"))
    if (identical(kind, "ordered")) {
      expect_equal(fitted(dist.profile), fitted(dist.dense),
                   tolerance = 1e-12, info = paste(kind, "distribution"))
      expect_equal(se(dist.profile), se(dist.dense),
                   tolerance = 1e-12, info = paste(kind, "distribution se"))
    }
  }

  run_case("unordered", 0.25, 2001L)
  run_case("ordered", 0.25, 2002L)
})

test_that("npudens categorical bandwidth search is unchanged by tree profile route", {
  old.tree <- getOption("np.tree")
  old.compress <- getOption("np.categorical.compress")
  old.messages <- getOption("np.messages")
  on.exit({
    options(np.tree = old.tree)
    options(np.categorical.compress = old.compress)
    options(np.messages = old.messages)
  }, add = TRUE)
  options(np.messages = FALSE)

  set.seed(20260608L)
  dat <- data.frame(
    x = ordered(rbinom(800L, 10L, 0.5))
  )

  options(np.tree = FALSE, np.categorical.compress = FALSE)
  dense <- npudensbw(~ x, data = dat, nmulti = 1)

  options(np.tree = FALSE, np.categorical.compress = TRUE)
  profile <- npudensbw(~ x, data = dat, nmulti = 1)

  expect_equal(profile$bw, dense$bw, tolerance = 1e-9)
  expect_equal(profile$fval, dense$fval, tolerance = 1e-10)
})

test_that("npudens categorical bootstrap profile route matches dense count algebra", {
  boot_fun <- getFromNamespace(".np_inid_boot_from_ksum_unconditional", "np")

  old.tree <- getOption("np.tree")
  old.compress <- getOption("np.categorical.compress")
  old.messages <- getOption("np.messages")
  on.exit({
    options(np.tree = old.tree)
    options(np.categorical.compress = old.compress)
    options(np.messages = old.messages)
  }, add = TRUE)
  options(np.messages = FALSE)

  set.seed(20260609L)
  dat <- data.frame(
    x = ordered(rbinom(300L, 10L, 0.5))
  )
  edat <- data.frame(x = ordered(levels(dat$x), levels = levels(dat$x)))
  bw <- npudensbw(~ x, data = dat, edat = edat,
                  bws = 0.25, bandwidth.compute = FALSE)
  counts <- replicate(5L, tabulate(sample.int(nrow(dat), nrow(dat), TRUE),
                                   nbins = nrow(dat)))

  options(np.tree = FALSE, np.categorical.compress = FALSE)
  dense <- boot_fun(
    xdat = dat,
    exdat = edat,
    bws = bw,
    B = ncol(counts),
    operator = "normal",
    counts = counts
  )

  options(np.tree = FALSE, np.categorical.compress = TRUE)
  profile <- boot_fun(
    xdat = dat,
    exdat = edat,
    bws = bw,
    B = ncol(counts),
    operator = "normal",
    counts = counts
  )

  expect_equal(profile$t0, dense$t0, tolerance = 1e-12)
  expect_equal(profile$t, dense$t, tolerance = 1e-12)
})

test_that("npudens formula route accepts numeric bws for categorical data", {
  old.tree <- getOption("np.tree")
  old.compress <- getOption("np.categorical.compress")
  old.messages <- getOption("np.messages")
  on.exit({
    options(np.tree = old.tree)
    options(np.categorical.compress = old.compress)
    options(np.messages = old.messages)
  }, add = TRUE)
  options(np.tree = FALSE, np.categorical.compress = TRUE, np.messages = FALSE)

  set.seed(20260610L)
  x <- ordered(rbinom(300L, 10L, 0.5))
  fit <- npudens(~ x, bws = 0.25)

  expect_s3_class(fit, "npdensity")
  expect_length(fitted(fit), length(x))
})

test_that("npcdens categorical ML bandwidth search uses exact profile route", {
  old.tree <- getOption("np.tree")
  old.compress <- getOption("np.categorical.compress")
  old.messages <- getOption("np.messages")
  on.exit({
    options(np.tree = old.tree)
    options(np.categorical.compress = old.compress)
    options(np.messages = old.messages)
  }, add = TRUE)
  options(np.messages = FALSE)

  set.seed(20260611L)
  n <- 700L
  x <- ordered(rbinom(n, 10L, 0.5))
  p <- plogis(as.numeric(as.character(x)) / 5 - 1)
  dat <- data.frame(
    y = ordered(rbinom(n, 10L, p)),
    x = x
  )

  options(np.tree = FALSE, np.categorical.compress = FALSE)
  dense <- npcdensbw(y ~ x, data = dat, bwmethod = "cv.ml", nmulti = 1)

  options(np.tree = FALSE, np.categorical.compress = TRUE)
  profile <- npcdensbw(y ~ x, data = dat, bwmethod = "cv.ml", nmulti = 1)

  expect_equal(profile$ybw, dense$ybw, tolerance = 1e-8)
  expect_equal(profile$xbw, dense$xbw, tolerance = 1e-8)
  expect_equal(profile$fval, dense$fval, tolerance = 1e-8)
})

test_that("npcdens categorical LS bandwidth search uses exact profile objective", {
  old.tree <- getOption("np.tree")
  old.compress <- getOption("np.categorical.compress")
  old.messages <- getOption("np.messages")
  on.exit({
    options(np.tree = old.tree)
    options(np.categorical.compress = old.compress)
    options(np.messages = old.messages)
  }, add = TRUE)
  options(np.messages = FALSE)

  set.seed(20260612L)
  n <- 450L
  x1 <- ordered(rbinom(n, 10L, 0.5))
  x2 <- factor(rbinom(n, 1L, 0.45))
  p <- plogis(as.numeric(as.character(x1)) / 5 + as.numeric(x2) / 3 - 1)
  dat <- data.frame(
    y = ordered(rbinom(n, 10L, p)),
    x1 = x1,
    x2 = x2
  )

  options(np.tree = FALSE, np.categorical.compress = FALSE)
  dense <- npcdensbw(y ~ x1 + x2, data = dat, bwmethod = "cv.ls", nmulti = 1)

  options(np.tree = FALSE, np.categorical.compress = TRUE)
  profile <- npcdensbw(y ~ x1 + x2, data = dat, bwmethod = "cv.ls", nmulti = 1)

  expect_equal(profile$ybw, dense$ybw, tolerance = 1e-5)
  expect_equal(profile$xbw, dense$xbw, tolerance = 1e-5)
  expect_equal(profile$fval, dense$fval, tolerance = 1e-8)
})

test_that("npcdist categorical bandwidth search uses exact profile objective", {
  old.tree <- getOption("np.tree")
  old.compress <- getOption("np.categorical.compress")
  old.messages <- getOption("np.messages")
  on.exit({
    options(np.tree = old.tree)
    options(np.categorical.compress = old.compress)
    options(np.messages = old.messages)
  }, add = TRUE)
  options(np.messages = FALSE)

  set.seed(20260613L)
  n <- 450L
  x1 <- ordered(rbinom(n, 10L, 0.5))
  x2 <- factor(rbinom(n, 1L, 0.45))
  p <- plogis(as.numeric(as.character(x1)) / 5 + as.numeric(x2) / 3 - 1)
  dat <- data.frame(
    y = ordered(rbinom(n, 10L, p)),
    x1 = x1,
    x2 = x2
  )

  options(np.tree = FALSE, np.categorical.compress = FALSE)
  dense <- npcdistbw(y ~ x1 + x2, data = dat, nmulti = 1)

  options(np.tree = FALSE, np.categorical.compress = TRUE)
  profile <- npcdistbw(y ~ x1 + x2, data = dat, nmulti = 1)

  expect_equal(profile$ybw, dense$ybw, tolerance = 1e-5)
  expect_equal(profile$xbw, dense$xbw, tolerance = 1e-5)
  expect_equal(profile$fval, dense$fval, tolerance = 1e-8)
})
