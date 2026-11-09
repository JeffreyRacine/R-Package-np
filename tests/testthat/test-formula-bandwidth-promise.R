test_that("formula selector handoff reuses a numeric bandwidth promise", {
  skip_on_cran()
  spawn_mpi_slaves(1L)
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  d <- data.frame(x = seq(.02, .98, length.out = 20))
  d$y <- sin(4 * d$x)
  hits <- 0L
  bandwidth <- function() { hits <<- hits + 1L; .3 * hits }
  fit <- npreg(y ~ x, data = d, bws = bandwidth(), bandwidth.compute = FALSE)
  expect_identical(hits, 1L)
  expect_identical(fit$bws$call[["..."]][["bws"]], quote(bandwidth()))
  expect_equal(as.numeric(fit$bws$bw), .3, tolerance = 0)
  ref <- npreg(y ~ x, data = d, bws = .3, bandwidth.compute = FALSE)
  expect_equal(fitted(fit), fitted(ref), tolerance = 1e-12)
})


 test_that("numeric bandwidth roles preserve forwarding and positional dots", {
   skip_on_cran()
  spawn_mpi_slaves(1L)
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
   d <- data.frame(x = seq(.02, .98, length.out = 20))
   d$y <- sin(4 * d$x)
   forward <- function(...) npreg(...)
   for (style in c("partial", "forwarded")) {
     hits <- 0L
     bw <- function() { hits <<- hits + 1L; .3 * hits }
     fit <- if (style == "partial")
       npreg(y ~ x, data = d, bw = bw(), bandwidth.compute = FALSE) else
       forward(y ~ x, data = d, bws = bw(), bandwidth.compute = FALSE)
     expect_identical(hits, 1L)
     expect_equal(as.numeric(fit$bws$bw), .3, tolerance = 0)
   }
   probe <- function(bws, xdat, ...) {
     .np_bw_formula_promises(
       quote(selector(y ~ x, d, NULL, na.omit, NULL, six, seven)), 2L,
       c(bws = "bws", xdat = "xdat"), environment(), parent.frame())
   }
   owned <- probe(y ~ x, d, NULL, na.omit, NULL, 6, 7)
   expect_identical(owned$call[[8L]], quote(seven))
   expect_false("bws" %in% names(owned$expressions))
 })
