test_that("conditional metadata validates physical categorical bandwidths", {
  withr::local_options(list(np.messages = FALSE, np.tree = FALSE))
  run <- function() {
    x <- data.frame(u = factor(rep(c("a", "b"), 8)),
                    o = ordered(rep(c("a", "b"), each = 8)))
    y <- factor(rep(c("a", "b"), 8))
    for (owner in c("npcdensbw", "npcdistbw")) {
      b <- getFromNamespace(owner, "npRmpi")(
        xdat = x, ydat = y, bws = c(1.25, 1.5, 3),
        bwscaling = TRUE, bandwidth.compute = FALSE,
        regtype = "lc", nomad = FALSE
      )
      expect_identical(b$ncatfac, 0.25)
      expect_equal(unname(b$bandwidth$x), c(.375, .75), tolerance = 0)
      expect_equal(unname(b$bandwidth$y), .3125, tolerance = 0)
      expect_equal(unname(b$sfactor$x), c(1.5, 3), tolerance = 0)
    }
  }
  getFromNamespace(".npRmpi_with_local_regression", "npRmpi")(run())
})

test_that("NN metadata retains counts when scaling is requested", {
  withr::local_options(list(np.messages = FALSE, np.tree = FALSE))
  run <- function() {
    i <- 1:16
    x <- data.frame(x = sin(i*sqrt(2)) + i/16)
    y <- cos(i*sqrt(3))
    for (owner in c("npregbw", "npudensbw", "npudistbw", "npcdensbw", "npcdistbw")) {
      args <- list(bws = 2, bwtype = "generalized_nn", bwscaling = TRUE,
                   bandwidth.compute = FALSE)
      if (owner %in% c("npudensbw", "npudistbw")) args$dat <- x
      else {
        args$xdat <- x; args$ydat <- y; args$regtype <- "lc"; args$nomad <- FALSE
        if (owner %in% c("npcdensbw", "npcdistbw")) args$bws <- c(2, 2)
      }
      b <- do.call(getFromNamespace(owner, "npRmpi"), args)
      expect_true(all(unlist(b$bandwidth) == 2))
      expect_true(all(unlist(b$sfactor) == 2))
    }
  }
  getFromNamespace(".npRmpi_with_local_regression", "npRmpi")(run())
})

test_that("categorical transforms preserve scale-factor units", {
  withr::local_options(list(np.messages = FALSE, np.tree = FALSE))
  run <- function() {
    x <- data.frame(o = ordered(rep(c("a", "b"), 8)))
    y <- sin((1:16)*sqrt(2)) + as.integer(x$o)/3
    b <- npregbw(xdat = x, ydat = y, bws = 3, bwscaling = TRUE,
      bandwidth.compute = FALSE, regtype = "lc", nomad = FALSE)
    core <- getFromNamespace(".npregbw_call_fixed_degree_core", "npRmpi")
    plain <- core(xdat = x, ydat = y, bws = b, eval.only = TRUE,
      transform.bounds = FALSE, invalid.penalty = "dbmax")
    mapped <- core(xdat = x, ydat = y, bws = b, eval.only = TRUE,
      transform.bounds = TRUE, invalid.penalty = "dbmax")
    expect_identical(mapped$bw, plain$bw)
    expect_identical(mapped$fval, plain$fval)
    expect_identical(mapped$num.feval, plain$num.feval)
  }
  getFromNamespace(".npRmpi_with_local_regression", "npRmpi")(run())
})

test_that("raw objective initialization preserves valid categorical scale factors", {
  withr::local_options(list(np.messages = FALSE, np.tree = FALSE))
  run <- function() {
    x <- data.frame(o = ordered(rep(c("a", "b"), 8)))
    y <- sin((1:16)*sqrt(2)) + as.integer(x$o)/3
    make <- function(scaling, value) npregbw(
      xdat = x, ydat = y, bws = value, bwscaling = scaling,
      bandwidth.compute = FALSE, regtype = "lc", nomad = FALSE
    )
    scaled <- make(TRUE, 3)
    physical <- make(FALSE, .75)
    evaluate <- getFromNamespace(".npregbw_eval_only", "npRmpi")
    s <- evaluate(xdat = x, ydat = y, bws = scaled, invalid.penalty = "dbmax")
    p <- evaluate(xdat = x, ydat = y, bws = physical, invalid.penalty = "dbmax")
    expect_true(is.finite(p$objective))
    expect_lt(p$objective, .Machine$double.xmax)
    expect_identical(s$objective, p$objective)
  }
  getFromNamespace(".npRmpi_with_local_regression", "npRmpi")(run())
})
