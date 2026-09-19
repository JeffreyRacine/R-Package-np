test_that("partial formula responses own retained rows and new omissions", {
  if (exists("spawn_mpi_slaves", mode = "function")) {
    spawn_mpi_slaves()
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }
  set.seed(8193)
  d <- data.frame(y = rnorm(24), x = runif(24), z = runif(24))
  rownames(d) <- paste0("subject-", 1:24)
  d$y[3] <- NA_real_
  visits <- 0L
  once <- function(x) { visits <<- visits + 1L; x^2 }
  b <- npregbw(y ~ once(x) + z, data = d, subset = c(24:1, 5L),
                bws = c(.4, .4), bandwidth.compute = FALSE, na.action = na.exclude)
  before <- b[[".np.formula.training"]]$frame
  frame <- b[[".np.formula.training"]]$frame
  n <- nrow(frame)
  yy <- seq_len(n)/n
  yy[4] <- NA_real_
  count <- visits
  seed <- .Random.seed
  actual <- npreg(b, tydat = yy)
  expect_identical(visits, count)
  expect_identical(.Random.seed, seed)
  expect_identical(b[[".np.formula.training"]]$frame, before)
  oracle <- npreg(b, txdat = frame[c("once(x)", "z")][-4, ], tydat = yy[-4])
  expect_equal(as.numeric(fitted(actual))[-4], as.numeric(fitted(oracle)), tolerance = 1e-10)
  expect_true(is.na(fitted(actual)[4]))
  expect_equal(as.integer(actual$omit), 4L)
  expect_equal(nrow(actual$bws[[".np.formula.training"]]$frame), n - 1L)
  expect_equal(as.numeric(fitted(npreg(actual$bws))), as.numeric(fitted(actual)))
  expect_error(npreg(b, tydat = yy[-1]), "retained training observation")
  expect_error(npreg(b, tydat = yy, data = d), "cannot be combined")
  expect_error(npreg(b, tydat = yy, na.action = na.fail), "missing values")
  clean <- npreg(b, tydat = replace(yy, 4, 0), na.action = na.omit)
  expect_length(fitted(clean), n)
  expect_null(clean$omit)
  # Replacing the evaluated response must not apply its transform again.
  t <- npregbw(log(y) ~ x, data = transform(d, y = exp(x)),
                bws = .4, bandwidth.compute = FALSE)
  rt <- npreg(t, tydat = seq_len(24)/24)
  expect_equal(fitted(rt), fitted(npreg(t, txdat = d["x"], tydat = seq_len(24)/24)))
})

test_that("formula response replacement is shared by the estimator adapters", {
  skip_on_cran()
  if (exists("spawn_mpi_slaves", mode = "function")) {
    spawn_mpi_slaves()
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }
  set.seed(93841)
  d <- data.frame(y = rnorm(28), x = runif(28), z = runif(28))
  yy <- d$y + d$x
  for (family in c("npreg", "npcdens", "npcdist", "npqreg", "npindex", "npplreg", "npscoef")) {
    ctor <- get(paste0(if (family == "npqreg") "npcdist" else family, "bw"))
    fit <- get(family)
    semi <- family %in% c("npplreg", "npscoef")
    cond <- family %in% c("npcdens", "npcdist", "npqreg")
    formula <- if (semi) y ~ x | z else y ~ x + z
    h <- if (family == "npplreg") matrix(.5, 2, 1) else if (family == "npscoef") .5 else
      if (family == "npindex") c(.5, 1, .4) else if (cond) rep(.5, 3) else c(.5, .5)
    b <- ctor(formula, data = d, bws = h, bandwidth.compute = FALSE)
    native <- list(txdat = if (semi) d["x"] else d[c("x", "z")], tydat = yy)
    if (semi) native$tzdat <- d["z"]
    oracle <- do.call(fit, c(list(bws = b, se = FALSE), native))
    actual <- fit(b, tydat = yy, se = FALSE)
    expect_equal(fitted(actual), fitted(oracle), tolerance = 1e-10, info = family)
    expect_equal(fitted(fit(actual$bws, se = FALSE)), fitted(actual), tolerance = 1e-10, info = family)
    nd <- transform(d[1:3, ], y = yy[1:3])
    ea <- if (semi) list(exdat = nd["x"], ezdat = nd["z"]) else
      list(exdat = nd[c("x", "z")])
    if (cond && family != "npqreg") ea$eydat <- nd["y"]
    expect_equal(fitted(fit(b, tydat = yy, newdata = nd, se = FALSE)),
                 fitted(do.call(fit, c(list(bws = b, se = FALSE), native, ea))),
                 tolerance = 1e-10, info = family)
    absent <- c(4L, 20L)
    changed <- fit(b, tydat = replace(yy, absent, NA_real_),
                   na.action = na.exclude, se = FALSE)
    retained <- lapply(native, function(x) if (is.data.frame(x))
      x[-absent, , drop = FALSE] else x[-absent])
    ref <- do.call(fit, c(list(bws = b, se = FALSE), retained))
    expect_equal(as.numeric(fitted(changed))[-absent], as.numeric(fitted(ref)),
                 tolerance = 1e-10, info = family)
    expect_true(all(is.na(fitted(changed)[absent])), info = family)
    expect_equal(fitted(fit(changed$bws, se = FALSE)), fitted(changed),
                 tolerance = 1e-10, info = family)
  }
  d$g <- factor(rep(c("a", "b"), 14))
  b <- npcdensbw(g ~ x, data = d, bws = c(.5, .2), bandwidth.compute = FALSE)
  gy <- factor(rep(c("b", "a"), 14), levels = levels(d$g))
  m <- npconmode(b, tydat = gy)
  expect_equal(fitted(m), fitted(npconmode(b, txdat = d["x"], tydat = gy)))
  expect_equal(fitted(npconmode(m$bws)), fitted(m))
  b <- npregbw(y ~ x, data = d, bws = .5, bandwidth.compute = FALSE)
  Y <- cbind(yy, yy^2)
  h <- npreghat(b, y = Y, output = "apply")
  ho <- npreghat(b, txdat = d["x"], y = Y, output = "apply")
  expect_equal(as.matrix(h), as.matrix(ho), ignore_attr = TRUE)
  Y[3, 2] <- NA_real_
  h <- npreghat(b, y = Y, output = "apply", na.action = na.exclude)
  ho <- npreghat(b, txdat = d[-3, "x", drop = FALSE], y = Y[-3, ], output = "apply")
  expect_equal(as.matrix(h), as.matrix(ho), ignore_attr = TRUE)
})
