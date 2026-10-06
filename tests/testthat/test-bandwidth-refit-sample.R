refit_sample_cases <- function() {
  set.seed(606)
  x <- data.frame(x = rnorm(120))
  z <- data.frame(z = rnorm(120))
  y <- .3 * x$x + sin(z$z) + rnorm(120, sd = .2)
  list(
    npregbw = list(xdat = x, ydat = y, bws = .7),
    npudensbw = list(dat = x, bws = .7),
    npudistbw = list(dat = x, bws = .7),
    npcdensbw = list(xdat = x, ydat = data.frame(y = y), bws = c(.7, .7)),
    npcdistbw = list(xdat = x, ydat = data.frame(y = y), bws = c(.7, .7)),
    npscoefbw = list(xdat = x, ydat = y, zdat = z, bws = .7),
    npplregbw = list(xdat = x, ydat = y, zdat = z, bws = matrix(.7, 2, 1)))
}

refit_sample_missing <- function(args) {
  roles <- intersect(names(args), c("dat", "xdat", "ydat", "zdat"))
  for (j in seq_along(roles)) {
    role <- roles[j]
    if (is.data.frame(args[[role]])) {
      args[[role]] <- args[[role]][seq_len(95), , drop = FALSE]
      args[[role]][3L * j, 1L] <- NA_real_
    } else {
      args[[role]] <- args[[role]][seq_len(95)]
      args[[role]][3L * j] <- NA_real_
    }
  }
  args
}

refit_sample_expect <- function(out, args) {
  roles <- intersect(names(args), c("dat", "xdat", "ydat", "zdat"))
  keep <- complete.cases(do.call(data.frame, args[roles]))
  expect_equal(out$nobs, sum(keep))
  expect_equal(out$nobs.omit, sum(!keep))
  expect_identical(as.integer(out$rows.omit), which(!keep))
  if (inherits(out, "plbandwidth"))
    expect_equal(vapply(out$bw, function(b) b$nobs, numeric(1)),
                 rep(sum(keep), length(out$bw)), ignore_attr = TRUE)
}

test_that("bandwidth refits use the current joint complete-case sample", {
  old.options <- options(np.messages = FALSE)
  on.exit(options(old.options), add = TRUE)
  cases <- refit_sample_cases()
  for (family in names(cases)) {
    args <- cases[[family]]
    original <- do.call(family, c(args, list(bandwidth.compute = FALSE)))
    frozen <- serialize(original, NULL)
    args <- refit_sample_missing(args)
    args$bws <- original
    out <- do.call(family, c(args, list(bandwidth.compute = FALSE)))
    refit_sample_expect(out, args)
    expect_identical(serialize(original, NULL), frozen)
  }
})

test_that("fixed-degree MADS preserves paired omission information", {
  skip_on_cran()
  old.options <- options(np.messages = FALSE)
  on.exit(options(old.options), add = TRUE)
  cases <- refit_sample_cases()
  for (family in c("npregbw", "npcdensbw", "npcdistbw")) {
    args <- refit_sample_missing(cases[[family]])
    out <- do.call(family, c(args, list(bwsolver = "mads", nmulti = 1L,
      nomad.opts = list(MAX_BB_EVAL = 30L), powell.remin = FALSE)))
    refit_sample_expect(out, args)
  }
})

test_that("automatic degree searches prepare one sample and retain its provenance", {
  skip_on_cran()
  old.options <- options(np.messages = FALSE)
  on.exit(options(old.options), add = TRUE)
  cases <- refit_sample_cases()
  for (family in c("npregbw", "npcdensbw", "npcdistbw", "npscoefbw", "npplregbw")) {
    args <- refit_sample_missing(cases[[family]])
    controls <- list(regtype = "lp", degree.select = "exhaustive",
      search.engine = "nomad", degree.min = 0L, degree.max = 1L,
      nmulti = 1L, powell.remin = FALSE, nomad.opts = list(MAX_BB_EVAL = 30L))
    if (family == "npscoefbw") controls$cv.iterate <- FALSE
    out <- do.call(family, c(args, controls))
    refit_sample_expect(out, args)
    roles <- intersect(names(args), c("xdat", "ydat", "zdat"))
    expect_identical(out[[".np.native.training", exact = TRUE]], args[roles])
  }
})
