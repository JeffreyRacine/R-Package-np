test_that("saved formula readers resolve data in its constructor owner", {
  if (!spawn_mpi_slaves()) skip("MPI test pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(1021L)
  d <- data.frame(y = rnorm(36), x = runif(36), z = rnorm(36))
  f <- y ~ x
  fs <- y ~ x | z
  make <- function(local.dat) npregbw(f, data = local.dat,
    bws = .6, bandwidth.compute = FALSE)
  make.expression <- function(local.dat) npregbw(f, data = local.dat[, c("x", "y")],
    bws = .6, bandwidth.compute = FALSE)
  make.sc <- function(local.dat) npscoefbw(fs, data = local.dat,
    bws = .6, bandwidth.compute = FALSE)
  canonical <- function(x) { attr(x, "call") <- NULL; x }
  test.fields <- c("In", "In.bootstrap", "P", "ixvar", "pivot.effective",
    "bootstrap.executed", "bootstrap.reason")
  for (make.bw in list(make, make.expression)) {
    bw <- make.bw(d)
    for (shadow in c(FALSE, TRUE)) {
      if (shadow) local.dat <- transform(d, x = x + 10, y = y * 4)
      h <- npreghat(bw)
      ref <- npreghat(bw, txdat = d["x"], y = d$y)
      expect_identical(canonical(h), canonical(ref))
      a <- npsigtest(bw, B = 9L, random.seed = 1022L)
      rng <- .Random.seed
      b <- npsigtest(bw, xdat = d["x"], ydat = d$y, B = 9L, random.seed = 1022L)
      expect_identical(a[test.fields], b[test.fields])
      expect_identical(.Random.seed, rng)
      override <- transform(d, x = x * 1.1, y = y + x)
      expect_identical(canonical(npreghat(bw, data = override)),
        canonical(npreghat(bw, txdat = override["x"], y = override$y)))
      a <- npsigtest(bw, data = override, B = 9L, random.seed = 1022L)
      rng <- .Random.seed
      b <- npsigtest(bw, xdat = override["x"], ydat = override$y,
        B = 9L, random.seed = 1022L)
      expect_identical(a[test.fields], b[test.fields])
      expect_identical(.Random.seed, rng)
      if (shadow) rm(local.dat)
    }
  }
  sb <- make.sc(d)
  local.dat <- transform(d, x = x + 10, y = y * 4)
  a <- plot(sb, output = "data", errors = "none", neval = 5L)
  b <- plot(sb, xdat = d["x"], ydat = d$y, zdat = d["z"],
    output = "data", errors = "none", neval = 5L)
  expect_identical(a, b)
  expect_error(npreghat(make(d), data = data.frame(wrong = 1:36)), "object 'y' not found")
  expect_length(fitted(npreg(make(d))), nrow(d))
})

test_that("saved readers retain subset, NA, transformed and aligned formula samples", {
  if (!spawn_mpi_slaves()) skip("MPI test pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(1023L)
  d <- data.frame(y = rnorm(40), x = runif(40), keep = rep(c(TRUE, TRUE, FALSE), length.out = 40))
  d$y[5L] <- NA_real_
  f <- y ~ x
  make <- function(local.dat) npregbw(f, data = local.dat, subset = keep,
    na.action = na.omit, bws = .6, bandwidth.compute = FALSE)
  bw <- make(d)
  mf <- model.frame(terms(bw), data = d, subset = keep, na.action = na.omit)
  new <- transform(d[1:8, ], x = x + .1)
  em <- model.frame(delete.response(terms(bw)), data = new)
  canonical <- function(x) { attr(x, "call") <- NULL; x }
  expect_identical(canonical(npreghat(bw, newdata = new)),
    canonical(npreghat(bw, txdat = mf[2L], y = mf[[1L]], exdat = em[1L])))
  fields <- c("In", "In.bootstrap", "P")
  a <- npsigtest(bw, B = 9L, random.seed = 1024L)
  rng <- .Random.seed
  b <- npsigtest(bw, xdat = mf[2L], ydat = mf[[1L]], B = 9L, random.seed = 1024L)
  expect_identical(a[fields], b[fields])
  expect_identical(.Random.seed, rng)
  expect_identical(a$rows.omit, as.vector(attr(mf, "na.action")))
  # Numeric transformed columns remain supported; matrix-valued poly
  # columns are rejected by the unchanged regression constructor.
  clean <- d[complete.cases(d), ]
  fp <- y ~ log(x)
  make.poly <- function(local.dat) npregbw(fp, data = local.dat,
    bws = .6, bandwidth.compute = FALSE)
  bp <- make.poly(clean)
  pm <- model.frame(terms(bp), data = clean)
  pe <- model.frame(delete.response(terms(bp)), data = new)
  expect_identical(canonical(npreghat(bp, newdata = new)),
    canonical(npreghat(bp, txdat = pm[2L], y = pm[[1L]], exdat = pe[1L])))
  y <- ts(rnorm(42), frequency = 4)
  ft <- y ~ lag(y, -1) + lag(y, -2)
  bt <- npregbw(ft, bws = c(.8, .9), bandwidth.compute = FALSE)
  tx <- data.frame(as.numeric(y)[2:41], as.numeric(y)[1:40])
  names(tx) <- c("lag(y, -1)", "lag(y, -2)")
  ty <- as.numeric(y)[3:42]
  expect_identical(canonical(npreghat(bt)), canonical(npreghat(bt, txdat = tx, y = ty)))
  a <- npsigtest(bt, B = 9L, random.seed = 1024L)
  rng <- .Random.seed
  b <- npsigtest(bt, xdat = tx, ydat = ty, B = 9L, random.seed = 1024L)
  expect_identical(a[fields], b[fields])
  expect_identical(.Random.seed, rng)
})

test_that("private MPI extraction respects saved data ownership", {
  if (!spawn_mpi_slaves()) skip("MPI test pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  set.seed(1025L)
  d <- data.frame(y = rnorm(36), x = runif(36))
  f <- y ~ x
  make <- function(local.dat) npregbw(f, data = local.dat,
    bws = .6, bandwidth.compute = FALSE)
  bw <- make(d)
  mf <- model.frame(f, data = d)
  expect_identical(.npRmpi_npsig_extract_xy_from_bws(bw),
    list(xdat = mf["x"], ydat = model.response(mf)))
})
