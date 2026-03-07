test_that("npsigtest basic functionality works", {
  set.seed(42)
  n <- 50 # Keep it small for speed
  x1 <- runif(n)
  x2 <- runif(n)
  y <- x1^2 + rnorm(n, sd=0.1) # x2 is irrelevant
  
  bw <- npregbw(y~x1+x2, bws=c(0.1, 0.5), bandwidth.compute=FALSE)
  
  # Significance test can be slow, use few boot replications
  sig <- npsigtest(bws=bw, boot.num=19)
  
  expect_s3_class(sig, "sigtest")
  expect_output(summary(sig))
})

test_that("npsigtest Type II hot-start helper honors nmulti contract", {
  helper <- getFromNamespace(".np_npsig_bootstrap_bw_reselect", "np")
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
                  bw.fun = stub_bw)
  helper(xdat = xdat,
         ydat = ydat,
         bws.seed = seed1,
         bootstrap.iter = 2L,
         extra.args = list(regtype = "ll"),
         bw.fun = stub_bw)

  expect_false("nmulti" %in% names(calls[[1L]]))
  expect_identical(calls[[1L]]$bws$seed_id, 0L)
  expect_identical(calls[[2L]]$nmulti, 1L)
  expect_identical(calls[[2L]]$bws$seed_id, 1L)
})

test_that("npsigtest Type II hot-start helper preserves explicit nmulti", {
  helper <- getFromNamespace(".np_npsig_bootstrap_bw_reselect", "np")
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
         bw.fun = stub_bw)
  helper(xdat = xdat,
         ydat = ydat,
         bws.seed = seed0,
         bootstrap.iter = 2L,
         extra.args = list(nmulti = 7L),
         bw.fun = stub_bw)

  expect_identical(calls[[1L]]$nmulti, 7L)
  expect_identical(calls[[2L]]$nmulti, 7L)
})

test_that("npsigtest formula interface path works", {
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

test_that("npsigtest Type II remains functional with hot-start reselection", {
  set.seed(17)
  n <- 30
  x1 <- runif(n)
  x2 <- runif(n)
  y <- x1 + rnorm(n, sd = 0.15)

  bw <- npregbw(y ~ x1 + x2, bws = c(0.2, 0.4), bandwidth.compute = FALSE)
  sig <- npsigtest(bws = bw, boot.num = 9, boot.type = "II")

  expect_s3_class(sig, "sigtest")
  expect_true(is.numeric(sig$P))
})

test_that("npsigtest rejects duplicate index entries", {
  set.seed(11)
  n <- 30
  x1 <- runif(n)
  x2 <- runif(n)
  y <- x1 + rnorm(n, sd = 0.1)

  bw <- npregbw(y ~ x1 + x2, bws = c(0.2, 0.4), bandwidth.compute = FALSE)

  expect_error(
    npsigtest(bws = bw, boot.num = 9, index = c(1, 1)),
    "repeated values"
  )
})
