test_that("native categorical coordinates replay public objectives exactly", {
  skip_if_not_installed("crs")
  old <- options(np.messages = FALSE, np.tree = FALSE,
                 np.categorical.compress = TRUE)
  on.exit(options(old), add = TRUE)

  set.seed(2324)
  n <- 53L
  unconditional <- data.frame(
    x = runif(n),
    u = factor(sample(letters[1:3], n, replace = TRUE)),
    o = ordered(sample(1:4, n, replace = TRUE), levels = 1:4)
  )
  x <- unconditional[c("x", "u", "o")]
  y <- data.frame(
    y = 0.35 * x$x + 0.14 * (x$u == "b") + rnorm(n, sd = 0.16)
  )
  nomad_options <- list(MAX_BB_EVAL = 14L)

  udens <- npudensbw(
    dat = unconditional, bwsolver = "mads", nmulti = 1L,
    nomad.opts = nomad_options
  )
  udens_replay <- np:::npudensbw.bandwidth(
    dat = unconditional, bws = udens, bandwidth.compute = TRUE,
    eval.only = TRUE, nmulti = 1L, invalid.penalty = "dbmax"
  )$fval[[1L]]

  udist_data <- unconditional[c("x", "o")]
  udist <- npudistbw(
    dat = udist_data, bwsolver = "mads", nmulti = 1L,
    bwscaling = FALSE, nomad.opts = nomad_options
  )
  udist_replay <- np:::npudistbw.dbandwidth(
    dat = udist_data, bws = udist, bandwidth.compute = TRUE,
    eval.only = TRUE, nmulti = 1L, invalid.penalty = "dbmax"
  )$fval[[1L]]

  cdens <- npcdensbw(
    xdat = x, ydat = y, bwsolver = "mads", nmulti = 1L,
    nomad.opts = nomad_options
  )
  cdens_replay <- np:::.npcdensbw_eval_only(
    x, y, cdens, invalid.penalty = "baseline"
  )$objective[[1L]]

  cdist <- npcdistbw(
    xdat = x, ydat = y, bwsolver = "mads", nmulti = 1L,
    bwscaling = FALSE, nomad.opts = nomad_options
  )
  cdist_replay <- np:::.npcdistbw_eval_only(
    x, y, bws = cdist, invalid.penalty = "baseline"
  )$objective[[1L]]

  expect_identical(as.numeric(udens$fval[[1L]]),
                   as.numeric(udens_replay))
  expect_identical(as.numeric(udist$fval[[1L]]),
                   as.numeric(udist_replay))
  expect_identical(as.numeric(cdens$fval[[1L]]),
                   as.numeric(cdens_replay))
  expect_identical(as.numeric(cdist$fval[[1L]]),
                   as.numeric(cdist_replay))
})
