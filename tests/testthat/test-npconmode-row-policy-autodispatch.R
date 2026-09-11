test_that("public mode fits transport row policy through actual autodispatch", {
  skip_on_cran()
  result <- npRmpi_run_isolated_contract(c(
    "library(npRmpi)",
    "npRmpi.init(nslaves = 1L, quiet = TRUE)",
    "options(np.messages = FALSE)",
    "x <- data.frame(x = seq(-1, 1, length.out = 48L))",
    "y <- data.frame(y = factor(rep(c('a', 'b'), 24L)))",
    "bw <- npcdensbw(xdat = x, ydat = y, bws = c(.15, .5),",
    "  bandwidth.compute = FALSE, cxkertype = 'epanechnikov', regtype = 'll')",
    "m <- npconmode(bws = bw, txdat = x, tydat = y, probabilities = TRUE)",
    "stopifnot(length(fitted(m)) == 48L, !anyNA(fitted(m)))",
    "ex <- data.frame(x = c(0, 4))",
    "state <- new.env(); state$n <- 0L",
    "e <- withCallingHandlers(npconmode(bws = bw, txdat = x, tydat = y,",
    "  exdat = ex, probabilities = TRUE, gradients = TRUE, se = TRUE),",
    "  warning = function(w) { state$n <- state$n + 1L; invokeRestart('muffleWarning') })",
    "stopifnot(!is.na(fitted(e)[1L]), is.na(fitted(e)[2L]), state$n == 1L)",
    "m2 <- npconmode(bws = bw, txdat = x, tydat = y)",
    "stopifnot(identical(fitted(m2), fitted(m)))",
    "cat('MODE_ROW_POLICY_AUTODISPATCH_OK\\n'); flush.console()"
  ), marker = "MODE_ROW_POLICY_AUTODISPATCH_OK", timeout = 25L)
  if (is.null(result)) skip("No installed npRmpi library for subprocess proof")
  expect_identical(result$status, 0L, info = paste(result$output, collapse = "\n"))
  expect_true(result$witnessed)
})
