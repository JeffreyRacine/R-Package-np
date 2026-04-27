#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(npRmpi))

route <- Sys.getenv("NP_RMPI_ATTACH_VALIDATE_ROUTE", "npreg")
route <- tolower(route)
valid.routes <- c("npreg", "npscoef", "npplreg", "npindex", "npcopula", "npconmode")
if (!(route %in% valid.routes)) {
  stop(
    sprintf(
      "NP_RMPI_ATTACH_VALIDATE_ROUTE must be one of: %s",
      paste(valid.routes, collapse = ", ")
    ),
    call. = FALSE
  )
}

is.master <- isTRUE(npRmpi.init(mode = "attach", quiet = TRUE))
options(np.messages = FALSE)

mark <- function(label) {
  cat(sprintf("ATTACH_ROUTE_STAGE %s\n", label), file = stderr())
  flush(stderr())
}

if (is.master) {
  suppressPackageStartupMessages(library(MASS))
  set.seed(42)
  n <- 80
  x <- runif(n)
  z <- runif(n)
  y <- sin(2 * pi * x) + 0.5 * z + rnorm(n, sd = 0.1)
  d <- data.frame(x = x, y = y)
  dsp <- data.frame(x = x, y = y, z = z)
  d2 <- data.frame(x = x, z = z)

  if (identical(route, "npreg")) {
    mark("npregbw")
    bw <- npregbw(y ~ x, data = d, regtype = "lc", bwmethod = "cv.ls", nmulti = 1)
    mark("npreg")
    fit <- npreg(bws = bw, data = d, gradients = FALSE)
    stopifnot(inherits(fit, "npregression"))
    cat("ATTACH_NPREG_ROUTE_OK\n")
  }

  if (identical(route, "npscoef")) {
    mark("npscoefbw")
    bw.sc <- npscoefbw(y ~ x | z, data = dsp, regtype = "lc", nmulti = 1)
    mark("npscoef")
    fit.sc <- npscoef(bws = bw.sc, data = dsp, gradients = FALSE)
    stopifnot(inherits(fit.sc, "smoothcoefficient"))
    cat("ATTACH_NPSCOEF_ROUTE_OK\n")
  }

  if (identical(route, "npplreg")) {
    mark("npplregbw")
    bw.pl <- npplregbw(y ~ x | z, data = dsp, regtype = "lc", nmulti = 1)
    mark("npplreg")
    fit.pl <- npplreg(bws = bw.pl, data = dsp, gradients = FALSE)
    stopifnot(inherits(fit.pl, "plregression"))
    cat("ATTACH_NPPLREG_ROUTE_OK\n")
  }

  if (identical(route, "npindex")) {
    mark("npindexbw")
    bw.si <- npindexbw(xdat = d2, ydat = y, regtype = "lc", nmulti = 1)
    mark("npindex")
    fit.si <- npindex(bws = bw.si, gradients = FALSE)
    stopifnot(inherits(fit.si, "singleindex"))
    cat("ATTACH_NPINDEX_ROUTE_OK\n")
  }

  if (identical(route, "npcopula")) {
    d.cop <- data.frame(x = x, y = y)
    u.cop <- data.frame(x = c(0.25, 0.5, 0.75), y = c(0.25, 0.5, 0.75))
    mark("npudistbw")
    bw.cop <- npudistbw(~x + y, data = d.cop)
    mark("npcopula")
    cop <- npcopula(bws = bw.cop, data = d.cop, u = u.cop, n.quasi.inv = 60)
    stopifnot(inherits(cop, "data.frame"))
    stopifnot(nrow(cop) == 9L)
    stopifnot(all(is.finite(cop$copula)))
    cat("ATTACH_NPCOPULA_ROUTE_OK\n")
  }

  if (identical(route, "npconmode")) {
    data(birthwt)
    bdat <- birthwt
    bdat$low <- factor(bdat$low)
    bdat$smoke <- factor(bdat$smoke)
    bdat$race <- factor(bdat$race)
    bdat$ht <- factor(bdat$ht)
    bdat$ui <- factor(bdat$ui)
    bdat$ftv <- ordered(bdat$ftv)
    mark("npcdensbw")
    bw.cm <- npcdensbw(low ~ smoke + race + ht + ui + ftv + age + lwt, data = bdat, nmulti = 1)
    mark("npconmode")
    fit.cm <- npconmode(bws = bw.cm)
    stopifnot(inherits(fit.cm, "conmode"))
    cat("ATTACH_NPCONMODE_ROUTE_OK\n")
  }

  cat("ATTACH_ROUTE_OK\n")
  npRmpi.quit(mode = "attach")
}
