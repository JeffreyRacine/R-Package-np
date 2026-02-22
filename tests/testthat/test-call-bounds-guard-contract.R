test_that(".Call wrappers tolerate empty continuous bounds", {
  set.seed(20260222)
  x <- rnorm(60)
  y <- x + rnorm(60)

  bw_reg <- npregbw(xdat = x, ydat = y, bws = 0.5, bandwidth.compute = FALSE, regtype = "lc")
  fit_reg <- npreg(bws = bw_reg)
  expect_s3_class(fit_reg, "npregression")

  bw_dens <- npudensbw(dat = data.frame(x = x), bws = 0.6, bandwidth.compute = FALSE)
  fit_dens <- npudens(tdat = data.frame(x = x), bws = bw_dens,
                      ckerlb = numeric(0), ckerub = numeric(0))
  expect_s3_class(fit_dens, "npdensity")

  bw_dist <- npudistbw(dat = data.frame(x = x), bws = 0.6, bandwidth.compute = FALSE)
  fit_dist <- npudist(tdat = data.frame(x = x), bws = bw_dist,
                      ckerlb = numeric(0), ckerub = numeric(0))
  expect_s3_class(fit_dist, "npdistribution")

  bw_cd <- npcdensbw(xdat = data.frame(x = x), ydat = data.frame(y = y),
                     bws = c(0.7, 0.5), bandwidth.compute = FALSE)
  fit_cd <- npcdens(
    txdat = data.frame(x = x), tydat = data.frame(y = y), bws = bw_cd,
    cxkerlb = numeric(0), cxkerub = numeric(0),
    cykerlb = numeric(0), cykerub = numeric(0)
  )
  expect_s3_class(fit_cd, "condensity")
})
