#!/usr/bin/env Rscript

# Repros for verified issues that were fixed in current working tree.
# Intended to be run against the updated package to confirm fixes.

suppressPackageStartupMessages(library(npRmpi))

cat("npRmpi version:", as.character(packageVersion("npRmpi")), "\n\n")

pool_ok <- TRUE
tryCatch({
  npRmpi.init(nslaves = 1, quiet = TRUE)
}, error = function(e) {
  pool_ok <<- FALSE
  cat("FAIL: npRmpi.init(nslaves=1) failed:", conditionMessage(e), "\n")
})

if (!pool_ok) {
  quit(save = "no", status = 1)
}

# Issue #4: cv.ls -> cv.ml segfault in npcdensbw with factor y
cat("#4: npcdensbw cv.ls -> cv.ml (factor y) ... ")
tryCatch({
  library(MASS)
  data(birthwt)
  birthwt$low <- factor(birthwt$low)
  nmulti <- 1
  bw1 <- npcdensbw(low ~ lwt, bwmethod = "cv.ls", data = birthwt, nmulti = nmulti)
  bw2 <- npcdensbw(low ~ lwt, bwmethod = "cv.ml", data = birthwt, nmulti = nmulti)
  cat("ok\n")
}, error = function(e) {
  cat("FAIL:", conditionMessage(e), "\n")
})

# Issue #6: bwtype propagation in npindexbw
cat("#6: npindexbw bwtype propagation ... ")
tryCatch({
  set.seed(123)
  n <- 120
  x1 <- runif(n, -1, 1)
  x2 <- runif(n, -1, 1)
  y <- ifelse(x1 + x2 + rnorm(n) > 0, 1, 0)
  bw_fixed <- npindexbw(y ~ x1 + x2, method = "kleinspady", gradients = TRUE, bwtype = "fixed")
  bw_adapt <- npindexbw(y ~ x1 + x2, method = "kleinspady", gradients = TRUE, bwtype = "adaptive_nn")
  bw_gen <- npindexbw(y ~ x1 + x2, method = "kleinspady", gradients = TRUE, bwtype = "generalized_nn")
  stopifnot(bw_fixed$type == "fixed")
  stopifnot(bw_adapt$type == "adaptive_nn")
  stopifnot(bw_gen$type == "generalized_nn")
  cat("ok\n")
}, error = function(e) {
  cat("FAIL:", conditionMessage(e), "\n")
})

# Issue #7: uocquantile factor subset mode
cat("#7: uocquantile factor subset ... ")
tryCatch({
  fruit.l <- factor(c("a","n","a","n","a","s"))
  m1 <- uocquantile(fruit.l, 0.5)
  m2 <- uocquantile(fruit.l[fruit.l != "a"], 0.5)
  stopifnot(as.character(m1) == "a")
  stopifnot(as.character(m2) == "n")
  cat("ok\n")
}, error = function(e) {
  cat("FAIL:", conditionMessage(e), "\n")
})

# Issue #13: npudens variance for categorical-only bws=0
cat("#13: npudens Avar for categorical-only ... ")
tryCatch({
  set.seed(42)
  n <- 100
  X <- sort(rbinom(n,2,0.4))
  p <- as.numeric(prop.table(table(X)))
  Avar.p <- p*(1-p)/n
  model <- npudens(tdat = factor(X), bws = 0, ukertype = "aitchisonaitken")
  p.hat <- unique(fitted(model))
  Avar.p.hat <- unique(se(model))^2
  stopifnot(all.equal(p, p.hat, tolerance = 1e-8))
  stopifnot(all.equal(Avar.p, Avar.p.hat, tolerance = 1e-8))
  cat("ok\n")
}, error = function(e) {
  cat("FAIL:", conditionMessage(e), "\n")
})

# Issue #5: npqreg tau validation
cat("#5: npqreg invalid tau ... ")
tryCatch({
  set.seed(1)
  x <- rnorm(30)
  y <- x + rnorm(30)
  ok <- FALSE
  tryCatch({
    npqreg(y ~ x, tau = 1.5)
  }, error = function(e) {
    ok <<- TRUE
  })
  stopifnot(ok)
  cat("ok\n")
}, error = function(e) {
  cat("FAIL:", conditionMessage(e), "\n")
})

# Issue #26: stored/dynamic formula in npplregbw
cat("#26: dynamic formula npplregbw ... ")
tryCatch({
  df <- data.frame(y=rnorm(10),x1=rnorm(10),x2=rnorm(10),x3=rnorm(10),x4=rnorm(10),x5=rnorm(10))
  f1 <- as.formula(paste0(paste("y ~ ", paste0("x", 1:4, collapse= " + ")), "|x5"))
  f2 <- y ~ x1 + x2 + x3 + x4 | x5
  npplregbw(formula = f1, data = df)
  npplregbw(formula = f2, data = df)
  cat("ok\n")
}, error = function(e) {
  cat("FAIL:", conditionMessage(e), "\n")
})

# Issue #51: npregiv exogenous covariates
cat("#51: npregiv with exogenous w ... ")
tryCatch({
  set.seed(123)
  y <- rnorm(60)
  x <- rnorm(60)
  z <- rnorm(60)
  w <- rnorm(60)
  npregiv(y, z, w, x)
  cat("ok\n")
}, error = function(e) {
  cat("FAIL:", conditionMessage(e), "\n")
})

# Session routing guard: npreg factor example from docs
cat("#reg: session npreg factor example ... ")
tryCatch({
  set.seed(42)
  n <- 250
  x <- runif(n)
  z1 <- rbinom(n, 1, 0.5)
  z2 <- rbinom(n, 1, 0.5)
  y <- cos(2 * pi * x) + z1 + rnorm(n, sd = 0.25)
  z1 <- factor(z1)
  z2 <- factor(z2)
  bw <- npregbw(y ~ x + z1 + z2, regtype = "lc", bwmethod = "cv.ls", nmulti = 1)
  fit <- npreg(bws = bw, gradients = FALSE)
  stopifnot(inherits(bw, "rbandwidth"), inherits(fit, "npregression"))
  cat("ok\n")
}, error = function(e) {
  cat("FAIL:", conditionMessage(e), "\n")
})

# Session routing guard: npcdens user-style example with plot()
cat("#npcdens: session user-style plot example ... ")
tryCatch({
  set.seed(42)
  n <- 1000
  x <- rnorm(n)
  y <- rnorm(n)
  F <- npcdens(y ~ x)
  summary(F$bws)
  png(tempfile(fileext = ".png"))
  on.exit(try(dev.off(), silent = TRUE), add = TRUE)
  plot(F)
  stopifnot(inherits(F, "condensity"))
  cat("ok\n")
}, error = function(e) {
  cat("FAIL:", conditionMessage(e), "\n")
})

# Session routing guard: npreghat parity vs npreg fit
cat("#npreghat: session fit-operator parity ... ")
tryCatch({
  set.seed(31415)
  n <- 120
  x <- runif(n)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.1)
  bw <- npregbw(y ~ x, regtype = "ll", bwmethod = "cv.ls", nmulti = 1)
  fit <- npreg(bws = bw, gradients = FALSE)
  H <- npreghat(bws = bw, txdat = data.frame(x = x))
  hy <- as.vector(H %*% y)
  stopifnot(inherits(H, "npreghat"))
  stopifnot(max(abs(hy - fitted(fit))) < 1e-6)
  cat("ok\n")
}, error = function(e) {
  cat("FAIL:", conditionMessage(e), "\n")
})

# Session routing guard: wild-hat bootstrap plot path
cat("#wildhat: session plot bootstrap path ... ")
tryCatch({
  set.seed(27182)
  n <- 80
  x <- runif(n)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.1)
  bw <- npregbw(y ~ x, bws = 0.2, bandwidth.compute = FALSE)
  old.chunk <- getOption("np.plot.wildhat.chunk.size")
  on.exit(options(np.plot.wildhat.chunk.size = old.chunk), add = TRUE)
  options(np.plot.wildhat.chunk.size = 5L)
  png(tempfile(fileext = ".png"))
  on.exit(try(dev.off(), silent = TRUE), add = TRUE)
  out <- suppressWarnings(plot(
    bw,
    xdat = data.frame(x = x),
    ydat = y,
    plot.behavior = "data",
    plot.errors.method = "bootstrap",
    plot.errors.boot.method = "wild-hat",
    plot.errors.boot.num = 9
  ))
  stopifnot(is.list(out), length(out) > 0)
  cat("ok\n")
}, error = function(e) {
  cat("FAIL:", conditionMessage(e), "\n")
})

cat("\nDone.\n")

# Issue #18: multi-dimensional instruments with Tikhonov should not error
cat("#18: npregiv Tikhonov with multi-dim w ... ")
tryCatch({
  set.seed(42)
  n <- 100
  v <- rnorm(n)
  eps <- rnorm(n, mean = 0, sd = 0.05)
  u <- -0.5 * v + eps
  w1 <- rnorm(n, mean = 0, sd = 1)
  w2 <- rnorm(n, mean = 0, sd = 1)
  w <- cbind(w1, w2)
  z <- 0.2 * w1 + v
  y <- z^2 + u
  npregiv(y = y, z = z, w = w, method = "Tikhonov")
  cat("ok\n")
}, error = function(e) {
  cat("FAIL:", conditionMessage(e), "\n")
})

# Deterministic teardown for non-interactive harness runs.
try(mpi.close.Rslaves(dellog = FALSE), silent = TRUE)
mpi.quit()
