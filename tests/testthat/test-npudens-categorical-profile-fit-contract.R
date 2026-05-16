test_that("npudens all-categorical tree profile fit matches dense fit", {
  old.tree <- getOption("np.tree")
  old.messages <- getOption("np.messages")
  on.exit({
    options(np.tree = old.tree)
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

    options(np.tree = FALSE)
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

    options(np.tree = TRUE)
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
  old.messages <- getOption("np.messages")
  on.exit({
    options(np.tree = old.tree)
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

    options(np.tree = FALSE)
    dens.bw.dense <- npudensbw(~ x, data = dat, edat = edat,
                               bws = bws, bandwidth.compute = FALSE)
    dens.dense <- npudens(bws = dens.bw.dense)
    if (identical(kind, "ordered")) {
      dist.bw.dense <- npudistbw(~ x, data = dat, edat = edat,
                                 bws = bws, bandwidth.compute = FALSE)
      dist.dense <- npudist(bws = dist.bw.dense)
    }

    options(np.tree = TRUE)
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
