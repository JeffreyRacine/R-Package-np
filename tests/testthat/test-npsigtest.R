test_that("npsigtest orchestrates locally without whole-call autodispatch", {
  fn <- getFromNamespace("npsigtest", "npRmpi")
  body.txt <- paste(deparse(body(fn), width.cutoff = 500L), collapse = " ")
  expect_false(grepl("\\.npRmpi_autodispatch_call\\(", body.txt))
  expect_false(grepl("\\.npRmpi_manual_distributed_call\\(", body.txt))
})

test_that("npsigtest Type II hot-start helper honors nmulti contract", {
  helper <- getFromNamespace(".npRmpi_npsig_bootstrap_bw_reselect", "npRmpi")
  xdat <- data.frame(x = c(0.1, 0.4, 0.8))
  ydat <- c(1, 2, 3)
  calls <- list()

  stub_bw <- function(...) {
    args <- list(...)
    calls[[length(calls) + 1L]] <<- args
    structure(list(bw = 0.25, seed_id = length(calls)), class = "rbandwidth")
  }

  seed0 <- structure(list(bw = 0.5, seed_id = 0L), class = "rbandwidth")
  seed1 <- helper(xdat = xdat,
                  ydat = ydat,
                  bws.seed = seed0,
                  bootstrap.iter = 1L,
                  extra.args = list(regtype = "ll"),
                  bw.fun = stub_bw,
                  localize = FALSE)
  helper(xdat = xdat,
         ydat = ydat,
         bws.seed = seed1,
         bootstrap.iter = 2L,
         extra.args = list(regtype = "ll"),
         bw.fun = stub_bw,
         localize = FALSE)

  expect_false("nmulti" %in% names(calls[[1L]]))
  expect_identical(calls[[1L]]$bws$seed_id, 0L)
  expect_identical(calls[[2L]]$nmulti, 1L)
  expect_identical(calls[[2L]]$bws$seed_id, 1L)
})

test_that("npsigtest Type II hot-start helper preserves explicit nmulti", {
  helper <- getFromNamespace(".npRmpi_npsig_bootstrap_bw_reselect", "npRmpi")
  xdat <- data.frame(x = c(0.1, 0.4, 0.8))
  ydat <- c(1, 2, 3)
  calls <- list()

  stub_bw <- function(...) {
    args <- list(...)
    calls[[length(calls) + 1L]] <<- args
    structure(list(bw = 0.25), class = "rbandwidth")
  }

  seed0 <- structure(list(bw = 0.5), class = "rbandwidth")
  helper(xdat = xdat,
         ydat = ydat,
         bws.seed = seed0,
         bootstrap.iter = 1L,
         extra.args = list(nmulti = 7L),
         bw.fun = stub_bw,
         localize = FALSE)
  helper(xdat = xdat,
         ydat = ydat,
         bws.seed = seed0,
         bootstrap.iter = 2L,
         extra.args = list(nmulti = 7L),
         bw.fun = stub_bw,
         localize = FALSE)

  expect_identical(calls[[1L]]$nmulti, 7L)
  expect_identical(calls[[2L]]$nmulti, 7L)
})

test_that("npsigtest gather helper preserves per-rank chunks when mpi.gather.Robj simplifies to a matrix", {
  gather.fun <- getFromNamespace(".npRmpi_npsig_gather_rank_chunks", "npRmpi")
  gathered <- matrix(1:9, nrow = 3L, ncol = 3L)

  chunks <- gather.fun(gathered = gathered, size = 3L)

  expect_type(chunks, "list")
  expect_identical(length(chunks), 3L)
  expect_identical(chunks[[1L]], 1:3)
  expect_identical(chunks[[2L]], 4:6)
  expect_identical(chunks[[3L]], 7:9)
})

test_that("npsigtest basic functionality works with autodispatch", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  options(npRmpi.autodispatch = TRUE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)

  set.seed(42)
  n <- 50 # Keep it small for speed
  x1 <- runif(n)
  x2 <- runif(n)
  y <- x1^2 + rnorm(n, sd=0.1) # x2 is irrelevant
  
  mydat <- data.frame(y, x1, x2)
  bw <- npregbw(y~x1+x2, data=mydat, bws=c(0.1, 0.5), bandwidth.compute=FALSE)
  
  # Significance test can be slow, use few boot replications
  sig <- npsigtest(bws=bw, boot.num=19)
  
  expect_s3_class(sig, "sigtest")
  expect_output(summary(sig))
})

test_that("npsigtest formula path works under autodispatch", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  options(npRmpi.autodispatch = TRUE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)

  set.seed(7)
  n <- 40
  x1 <- runif(n)
  x2 <- runif(n)
  y <- x1 + rnorm(n, sd = 0.1)
  mydat <- data.frame(y, x1, x2)

  sig <- npsigtest(y ~ x1 + x2,
                   data = mydat,
                   boot.num = 9)

  expect_s3_class(sig, "sigtest")
  expect_true(is.numeric(sig$P))
})

test_that("npsigtest npregression path works under autodispatch", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  options(npRmpi.autodispatch = TRUE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)

  set.seed(42)
  n <- 80
  z <- factor(rbinom(n, 1, .5))
  x1 <- rnorm(n)
  x2 <- runif(n, -2, 2)
  y <- x1 + x2 + rnorm(n, sd = 0.2)
  mydat <- data.frame(z, x1, x2, y)

  model <- npreg(y ~ z + x1 + x2,
                 regtype = "ll",
                 bwmethod = "cv.aic",
                 data = mydat)

  sig <- npsigtest(model, boot.num = 9)

  expect_s3_class(sig, "sigtest")
  expect_true(is.numeric(sig$P))
})

test_that("npsigtest rejects duplicate index entries under autodispatch", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  options(npRmpi.autodispatch = TRUE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)

  set.seed(11)
  n <- 30
  x1 <- runif(n)
  x2 <- runif(n)
  y <- x1 + rnorm(n, sd = 0.1)
  mydat <- data.frame(y, x1, x2)

  bw <- npregbw(y ~ x1 + x2, data = mydat, bws = c(0.2, 0.4), bandwidth.compute = FALSE)

  expect_error(
    npsigtest(bws = bw, boot.num = 9, index = c(1, 1)),
    "repeated values"
  )
})

test_that("npsigtest local regression wrapper is safe under an active slave pool", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  options(npRmpi.autodispatch = TRUE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)

  set.seed(21)
  n <- 35
  x <- runif(n)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.1)
  bw <- npregbw(xdat = data.frame(x = x), ydat = y, bws = 0.25, bandwidth.compute = FALSE, regtype = "ll")

  local.fun <- getFromNamespace(".npRmpi_npsig_npreg_local", "npRmpi")
  fit <- local.fun(txdat = data.frame(x = x), tydat = y, bws = bw, gradients = TRUE)

  expect_s3_class(fit, "regression")
  expect_identical(length(fit$mean), n)
  expect_identical(ncol(fit$grad), 1L)
})

test_that("npsigtest boot.type II works under autodispatch", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  options(npRmpi.autodispatch = TRUE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)

  set.seed(31)
  n <- 30
  x1 <- runif(n)
  x2 <- runif(n)
  y <- x1 + rnorm(n, sd = 0.1)
  mydat <- data.frame(y, x1, x2)
  bw <- npregbw(y ~ x1 + x2, data = mydat, bws = c(0.2, 0.4), bandwidth.compute = FALSE)

  sig <- npsigtest(bws = bw, boot.num = 9, boot.type = "II", joint = TRUE, index = 1)

  expect_s3_class(sig, "sigtest")
  expect_identical(length(sig$P), 1L)
})
