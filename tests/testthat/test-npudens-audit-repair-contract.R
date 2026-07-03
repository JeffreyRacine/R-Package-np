with_npudens_audit_runtime <- function(code) {
  code <- substitute(code)
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  if (exists("spawn_mpi_slaves", mode = "function")) {
    skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }

  eval(code, envir = parent.frame())
}

expect_omit <- function(x, value) {
  expect_identical(as.integer(x), as.integer(value))
}

cat_se_oracle <- function(dat) {
  bw <- npudensbw(dat = dat, bws = rep(0, ncol(dat)),
                  bandwidth.compute = FALSE)
  fit <- npudens(bws = bw, tdat = dat)
  p <- pmin(pmax(fit$dens, 0), 1)
  sqrt(p * (1 - p) / nrow(dat))
}

test_that("npudens formula reentry honors explicit data", {
  with_npudens_audit_runtime({
    d1 <- data.frame(x = seq(-1, 1, length.out = 8L))
    d2 <- data.frame(x = seq(-2, 2, length.out = 13L))

    bw <- npudensbw(~ x, data = d1, bws = 0.4,
                    bandwidth.compute = FALSE)
    fit <- npudens(bws = bw, data = d2)

    expect_identical(fit$nobs, nrow(d2))
    expect_identical(fit$ntrain, nrow(d2))
    expect_identical(length(fit$dens), nrow(d2))
  })
})

test_that("npudens keeps train and evaluation omitted rows separately", {
  with_npudens_audit_runtime({
    train <- data.frame(x = c(0, NA, 1, 2, NA, 3))
    eval.clean <- data.frame(x = c(-0.5, 0.5, 1.5))
    eval.dirty <- data.frame(x = c(-0.5, NA, 1.5, NA))
    bw <- npudensbw(dat = data.frame(x = c(0, 1, 2, 3)),
                    bws = 0.5, bandwidth.compute = FALSE)

    clean.fit <- npudens(bws = bw, tdat = train, edat = eval.clean)
    dirty.fit <- npudens(bws = bw, tdat = train, edat = eval.dirty)

    expect_omit(clean.fit$train.rows.omit, c(2L, 5L))
    expect_identical(clean.fit$ntrain.omit, 2L)
    expect_identical(clean.fit$eval.rows.omit, NA)
    expect_identical(clean.fit$neval.omit, 0)
    expect_identical(clean.fit$rows.omit, NA)
    expect_identical(clean.fit$nobs.omit, 0)

    expect_omit(dirty.fit$train.rows.omit, c(2L, 5L))
    expect_omit(dirty.fit$eval.rows.omit, c(2L, 4L))
    expect_omit(dirty.fit$rows.omit, c(2L, 4L))
    expect_identical(dirty.fit$ntrain.omit, 2L)
    expect_identical(dirty.fit$neval.omit, 2L)
    expect_identical(dirty.fit$nobs.omit, 2L)
  })
})

test_that("npudens zero-row inputs fail before low-level shape errors", {
  with_npudens_audit_runtime({
    bw <- npudensbw(dat = data.frame(x = c(0, 1, 2, 3)),
                    bws = 0.5, bandwidth.compute = FALSE)

    expect_error(
      npudens(bws = bw, tdat = data.frame(x = c(NA_real_, NA_real_))),
      "no rows without NAs"
    )
    expect_error(
      npudens(bws = bw, tdat = data.frame(x = c(0, 1, 2, 3)),
              edat = data.frame(x = c(NA_real_, NA_real_))),
      "no rows without NAs"
    )
    expect_error(
      npudensbw(dat = data.frame(x = c(NA_real_, NA_real_)),
                nmulti = 1L, itmax = 1L),
      "no rows without NAs"
    )
  })
})

test_that("npudens categorical exact-zero standard errors use binomial proportion scale", {
  with_npudens_audit_runtime({
    unordered <- data.frame(u = factor(c("a", "a", "b", "b", "b", "c")))
    ordered <- data.frame(o = ordered(c("low", "low", "mid", "mid", "high", "high"),
                                      levels = c("low", "mid", "high")))
    mixed <- data.frame(
      u = factor(c("a", "a", "a", "b", "b", "b")),
      o = ordered(c("low", "low", "mid", "mid", "mid", "high"),
                  levels = c("low", "mid", "high"))
    )

    for (dat in list(unordered, ordered, mixed)) {
      bw <- npudensbw(dat = dat, bws = rep(0, ncol(dat)),
                      bandwidth.compute = FALSE)
      fit <- npudens(bws = bw, tdat = dat)
      expect_equal(fit$derr, cat_se_oracle(dat), tolerance = 1e-14)
    }
  })
})

test_that("npudensbw validates method metadata without rejecting normal-reference", {
  with_npudens_audit_runtime({
    x <- data.frame(x = seq(0, 1, length.out = 8L))
    bw <- npudensbw(dat = x, bws = 0.3, bandwidth.compute = FALSE)
    bad <- bw
    bad$method <- "bogus"

    expect_error(
      npudensbw(dat = x, bws = bad, nmulti = 1L, itmax = 1L),
      "does not support bwmethod"
    )

    nr <- npudensbw(dat = data.frame(
      x = seq(0, 1, length.out = 8L),
      u = factor(rep(c("a", "b"), 4L)),
      o = ordered(rep(c("low", "high"), 4L), levels = c("low", "high"))
    ), bwmethod = "normal-reference")
    expect_identical(nr$method, "normal-reference")
    cat.mask <- rep(FALSE, length(nr$bw))
    if (isTRUE(nr$nuno > 0L))
      cat.mask <- cat.mask | nr$iuno
    if (isTRUE(nr$nord > 0L))
      cat.mask <- cat.mask | nr$iord
    expect_equal(nr$bw[cat.mask], c(0, 0), tolerance = 0)
  })
})

test_that("conditional bandwidth ordered liracine mappings remain normalized", {
  ns <- environment(npcdensbw)
  candidates <- ls(ns, all.names = TRUE)
  candidates <- candidates[grepl("npcdensbw|npcdistbw", candidates)]
  bodies <- vapply(candidates, function(name) {
    obj <- get(name, envir = ns)
    if (!is.function(obj))
      return("")
    paste(deparse(body(obj)), collapse = "\n")
  }, character(1L))

  bad <- names(bodies)[grepl("liracine = OKER_LR", bodies, fixed = TRUE)]
  expect_identical(bad, character(0))
  expect_true(any(grepl("liracine = OKER_NLR", bodies, fixed = TRUE)))
})
