library(np)

set.seed(70619)
n <- 30L
xdat <- data.frame(x = runif(n))
ydat <- rnorm(n)

bw <- npcdistbw(
  xdat = xdat,
  ydat = data.frame(y = ydat),
  bws = c(0.35, 0.35),
  bandwidth.compute = FALSE
)
fit <- npqreg(bws = bw, txdat = xdat, tydat = ydat, tau = 0.5)

for (boot in c("inid", "fixed", "geom", "wild")) {
  msg <- tryCatch({
    plot(
      fit,
      output = "data",
      perspective = FALSE,
      errors = "bootstrap",
      bootstrap = boot,
      center = "bias-corrected",
      B = 5L,
      neval = 4L
    )
    NA_character_
  }, error = conditionMessage)

  if (is.na(msg) ||
      !grepl("conditional quantile \\(npqreg\\) bootstrap plots", msg) ||
      !grepl("use center=\"estimate\"", msg, fixed = TRUE) ||
      !grepl(sprintf("bootstrap=\"%s\"", boot), msg, fixed = TRUE)) {
    stop(sprintf(
      "npqreg bias-center fail-closed sentinel failed for bootstrap=%s: %s",
      boot,
      msg
    ), call. = FALSE)
  }
}

cat("NPQREG_BIAS_CENTER_FAILCLOSED_OK\n")
