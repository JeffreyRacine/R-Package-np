#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(npRmpi))

iters <- as.integer(Sys.getenv("NP_NATIVE_STRESS_ITERS", "12"))
if (!is.finite(iters) || iters < 1L)
  stop("NP_NATIVE_STRESS_ITERS must be a positive integer")

options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = FALSE)

npRmpi.init(nslaves = 1, quiet = TRUE)
# In one-shot harness runs, process teardown is more reliable than explicit
# quit in this environment (explicit quit can hang after successful work).

set.seed(20260223)
n <- 120L
x <- runif(n)
y <- rnorm(n)

xdat <- data.frame(x = x)
ydat <- data.frame(y = y)

bw_reg <- npregbw(xdat = xdat, ydat = y, bws = 0.35, bandwidth.compute = FALSE, regtype = "lc")
bw_ud <- npudensbw(dat = xdat, bws = 0.35, bandwidth.compute = FALSE,
                   ckerlb = numeric(0), ckerub = numeric(0))
bw_uf <- npudistbw(dat = xdat, bws = 0.35, bandwidth.compute = FALSE,
                   ckerlb = numeric(0), ckerub = numeric(0))
bw_cd <- npcdensbw(xdat = xdat, ydat = ydat, bws = c(0.35, 0.35), bandwidth.compute = FALSE,
                   cxkerlb = numeric(0), cxkerub = numeric(0),
                   cykerlb = numeric(0), cykerub = numeric(0))
bw_cf <- npcdistbw(xdat = xdat, ydat = ydat, bws = c(0.35, 0.35), bandwidth.compute = FALSE,
                   cxkerlb = numeric(0), cxkerub = numeric(0),
                   cykerlb = numeric(0), cykerub = numeric(0))

for (i in seq_len(iters)) {
  fit_r <- npreg(bws = bw_reg, txdat = xdat, tydat = y, exdat = xdat, gradients = FALSE)
  fit_ud <- npudens(bws = bw_ud, tdat = xdat)
  fit_uf <- npudist(bws = bw_uf, tdat = xdat)
  fit_cd <- npcdens(bws = bw_cd, txdat = xdat, tydat = ydat)
  fit_cf <- npcdist(bws = bw_cf, txdat = xdat, tydat = ydat)
  fit_ks <- npksum(txdat = xdat, tydat = y, exdat = xdat, bws = 0.35,
                   ckerlb = numeric(0), ckerub = numeric(0))

  stopifnot(
    inherits(fit_r, "npregression"),
    inherits(fit_ud, "npdensity"),
    inherits(fit_uf, "npdistribution"),
    inherits(fit_cd, "condensity"),
    inherits(fit_cf, "condistribution"),
    inherits(fit_ks, "npkernelsum")
  )

  if (i %% 4L == 0L)
    cat("iter", i, "ok\n")
}

cat("NATIVE_BRIDGE_STRESS_NPRMPI_OK iters=", iters, "\n", sep = "")
flush.console()
quit(save = "no", status = 0, runLast = FALSE)
