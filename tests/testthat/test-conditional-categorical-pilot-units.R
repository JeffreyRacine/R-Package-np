test_that("conditional fallback categorical pilots preserve physical units", {
  package <- if ("npRmpi" %in% loadedNamespaces()) "npRmpi" else "np"
  withr::local_options(list(np.messages = FALSE, np.tree = FALSE))
  run <- function() {
    make <- getFromNamespace("npcdensbw", package)
    evaluate <- getFromNamespace(".npcdensbw_eval_only", package)
    i <- seq_len(16L)
    u <- factor(c(1L, rep(2L, 7L), rep(3L, 8L)))
    fixtures <- list(
      unordered = list(x = data.frame(u = u), y = data.frame(y = sin(i*sqrt(2))),
                       start = c(1, 0), pilot = c(2, 1/3)),
      ordered = list(x = data.frame(o = ordered(u)),
                     y = data.frame(y = sin(i*sqrt(2))),
                     start = c(1, 0), pilot = c(2, .5)),
      split = list(x = data.frame(u = u, o = ordered(u)),
                   y = data.frame(y = sin(i*sqrt(2)), u = u, o = ordered(u)),
                   start = c(1, 0, 0, 0, 0), pilot = c(2, 1/3, .5, 1/3, .5)))
    for (label in names(fixtures)) {
      f <- fixtures[[label]]
      construct <- function(values, scaling = FALSE)
        make(xdat = f$x, ydat = f$y, bws = values, bandwidth.compute = FALSE,
             bwscaling = scaling, regtype = "lc", nomad = FALSE)
      score <- function(b, mode)
        evaluate(xdat = f$x, ydat = f$y, bws = b, invalid.penalty = mode)$objective
      physical <- construct(f$start)
      scaled <- construct(c(physical$sfactor$y, physical$sfactor$x), TRUE)
      pilot <- construct(f$pilot)
      raw <- score(pilot, "dbmax")
      expected <- raw - (abs(raw) + 1) * 10
      expect_true(is.finite(raw) && abs(raw) < .Machine$double.xmax, info = label)
      expect_identical(score(physical, "dbmax"), -.Machine$double.xmax)
      expect_identical(score(scaled, "dbmax"), -.Machine$double.xmax)
      expect_equal(score(physical, "baseline"), expected, tolerance = 1e-13)
      expect_equal(score(scaled, "baseline"), expected, tolerance = 1e-13)
      scaled.pilot <- construct(c(pilot$sfactor$y, pilot$sfactor$x), TRUE)
      expect_equal(score(pilot, "baseline"), score(scaled.pilot, "baseline"),
                   tolerance = 1e-13)
      if (label == "unordered") {
        weights <- ifelse(outer(f$x$u, f$x$u, "=="), 2/3, 1/6)
        diag(weights) <- 0
        ky <- dnorm(outer(f$y$y, f$y$y, "-")/2)/2
        reference <- sum(log(rowSums(weights * ky)/rowSums(weights)))
        expect_equal(raw, reference, tolerance = 1e-13)
      }
    }
  }
  if (package == "npRmpi")
    getFromNamespace(".npRmpi_with_local_regression", package)(run())
  else run()
})

test_that("adaptive conditional fallback pilots share the categorical units", {
  package <- if ("npRmpi" %in% loadedNamespaces()) "npRmpi" else "np"
  withr::local_options(list(np.messages = FALSE, np.tree = FALSE,
                            np.extendednn = TRUE))
  run <- function() {
    i <- seq_len(16L)
    x <- data.frame(z = rep(1:8, each = 2),
                    u = factor(c(1L, rep(2L, 7L), rep(3L, 8L))))
    y <- sin(i*sqrt(2))
    for (owner in c("npcdensbw", "npcdistbw")) {
      make <- getFromNamespace(owner, package)
      evaluate <- getFromNamespace(paste0(".", owner, "_eval_only"), package)
      construct <- function(values, scaling = FALSE)
        make(xdat = x, ydat = y, bws = values, bandwidth.compute = FALSE,
             bwscaling = scaling, bwtype = "adaptive_nn",
             regtype = "lc", nomad = FALSE)
      score <- function(b, mode)
        evaluate(xdat = x, ydat = y, bws = b, invalid.penalty = mode)$objective
      physical <- construct(c(1, 1, 0))
      scaled <- construct(c(physical$sfactor$y, physical$sfactor$x), TRUE)
      raw <- score(construct(c(2, 2, 1/3)), "dbmax")
      sign <- if (owner == "npcdensbw") -1 else 1
      expected <- raw + sign*(abs(raw) + 1)*10
      expect_identical(score(physical, "dbmax"), sign*.Machine$double.xmax)
      expect_identical(score(scaled, "dbmax"), sign*.Machine$double.xmax)
      expect_equal(score(physical, "baseline"), expected, tolerance = 1e-13)
      expect_equal(score(scaled, "baseline"), expected, tolerance = 1e-13)
    }
  }
  if (package == "npRmpi")
    getFromNamespace(".npRmpi_with_local_regression", package)(run())
  else run()
})
