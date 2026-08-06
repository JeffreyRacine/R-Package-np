test_that("resident npcdens cache keys bandwidths and polynomial degree", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260805)
  prior_n <- 32L
  prior_x <- data.frame(
    x1 = runif(prior_n),
    x2 = runif(prior_n)
  )
  prior_y <- data.frame(y = factor(ifelse(
    prior_x$x1 + prior_x$x2 + rnorm(prior_n, sd = 0.15) > 1,
    "yes",
    "no"
  )))
  invisible(np::npcdensbw(
    xdat = prior_x,
    ydat = prior_y,
    regtype = "lc",
    bwmethod = "cv.ml",
    bwtype = "fixed",
    bwsolver = "powell",
    nmulti = 1L
  ))

  set.seed(271828)
  n <- 70L
  x <- sort(runif(n))
  y <- sin(2 * pi * x) + 0.3 * x + rnorm(n, sd = 0.18)
  xdat <- data.frame(x = x)
  ydat <- data.frame(y = y)
  bw <- np::npcdensbw(
    xdat = xdat,
    ydat = ydat,
    regtype = "lp",
    degree = 0L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    bws = c(0.24, 0.31),
    bandwidth.compute = FALSE
  )
  prep <- np:::.npcdensbw_prepared_prepare_args(
    xdat = xdat,
    ydat = ydat,
    bws = bw,
    start.bw = NULL,
    invalid.penalty = "baseline",
    penalty.multiplier = 10,
    degree.search = TRUE
  )
  points <- data.frame(
    xbw = c(0.24, 0.24, 0.41, 0.24),
    ybw = c(0.31, 0.31, 0.37, 0.31),
    degree = c(0L, 2L, 2L, 0L)
  )

  evaluate <- function(cache) {
    options(np.objective.cache = cache)
    do.call(np:::npPreparedObjectivePrepareConditionalDensity, list(
      c.uno = prep$c.uno, c.ord = prep$c.ord, c.con = prep$c.con,
      u.uno = prep$u.uno, u.ord = prep$u.ord, u.con = prep$u.con,
      mysd = prep$mysd, myopti = prep$myopti, myoptd = prep$myoptd,
      rbw = prep$rbw, penalty.mode = prep$penalty_mode,
      penalty.multiplier = prep$penalty_multiplier, degree = prep$degree,
      bernstein = prep$bernstein, basis = prep$basis, regtype = prep$regtype,
      cxkerlb = prep$cxkerlb, cxkerub = prep$cxkerub,
      cykerlb = prep$cykerlb, cykerub = prep$cykerub
    ))
    on.exit(np:::npPreparedObjectiveDestroyConditionalDensity(), add = TRUE)
    vapply(seq_len(nrow(points)), function(i) {
      np:::npPreparedObjectiveEvalConditionalDensity(
        bw = c(points$xbw[i], points$ybw[i]),
        degree = points$degree[i]
      )[1L]
    }, numeric(1L))
  }

  cached <- evaluate(TRUE)
  uncached <- evaluate(FALSE)
  direct <- vapply(seq_len(nrow(points)), function(i) {
    candidate <- bw
    candidate$regtype <- "lp"
    candidate$pregtype <- "Local-Polynomial"
    candidate$regtype.engine <- "lp"
    candidate$degree <- as.integer(points$degree[i])
    candidate$degree.engine <- as.integer(points$degree[i])
    candidate$xbw[candidate$ixcon] <- points$xbw[i]
    candidate$ybw[candidate$iycon] <- points$ybw[i]
    np:::.npcdensbw_eval_only(
      xdat = xdat,
      ydat = ydat,
      bws = candidate
    )$objective
  }, numeric(1L))

  expect_equal(cached, uncached, tolerance = 1e-13)
  expect_equal(cached, direct, tolerance = 1e-13)
  expect_false(isTRUE(all.equal(cached[1L], cached[2L], tolerance = 1e-8)))
  expect_identical(cached[1L], cached[4L])
})

test_that("prepared npcdens context bounds degree writes by allocated key length", {
  src <- readLines(np_test_source_path("src", "np.c"), warn = FALSE)
  text <- paste(src, collapse = "\n")
  expect_match(
    text,
    "np_conditional_density_prepared_context\\.degree_key_len = degree_key_len"
  )
  expect_match(
    text,
    "i < np_conditional_density_prepared_context\\.degree_key_len"
  )
  expect_match(
    text,
    "np_conditional_density_prepared_context.num_all_var \\+ i \\+ 1"
  )
})
