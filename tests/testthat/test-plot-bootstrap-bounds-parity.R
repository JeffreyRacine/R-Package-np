test_that("bootstrap quantile bounds preserve legacy numerical results", {
  skip_if_not_installed("npRmpi")

  new_bounds <- getFromNamespace("compute.bootstrap.quantile.bounds", "npRmpi")
  new_scs <- getFromNamespace("np.plot.SCSrank", "npRmpi")

  old_scs <- function(x, conf.level = 0.95, alternative = "two.sided") {
    alternative <- match.arg(alternative, choices = c("two.sided", "less", "greater"))
    DataMatrix <- x
    N <- nrow(DataMatrix)
    k <- round(conf.level * N, 0)
    RankDat <- apply(DataMatrix, 2, rank)

    switch(alternative,
      "two.sided" = {
        W1 <- apply(RankDat, 1, max)
        W2 <- N + 1 - apply(RankDat, 1, min)
        Wmat <- cbind(W1, W2)
        w <- apply(Wmat, 1, max)
        tstar <- round(sort(w)[k], 0)
        SCI <- function(x) {
          sortx <- sort(x)
          cbind(sortx[N + 1 - tstar], sortx[tstar])
        }
        SCS <- t(apply(DataMatrix, 2, SCI))
      },
      "less" = {
        W1 <- apply(RankDat, 1, max)
        tstar <- round(sort(W1)[k], 0)
        SCI <- function(x) {
          sortx <- sort(x)
          cbind(-Inf, sortx[tstar])
        }
        SCS <- t(apply(DataMatrix, 2, SCI))
      },
      "greater" = {
        W2 <- N + 1 - apply(RankDat, 1, min)
        tstar <- round(sort(W2)[k], 0)
        SCI <- function(x) {
          sortx <- sort(x)
          cbind(sortx[N + 1 - tstar], Inf)
        }
        SCS <- t(apply(DataMatrix, 2, SCI))
      }
    )

    colnames(SCS) <- c("lower", "upper")
    attr(SCS, which = "k") <- k
    attr(SCS, which = "N") <- N
    list(conf.int = SCS, conf.level = conf.level, alternative = alternative)
  }

  old_bounds <- function(boot.t, alpha, band.type) {
    B <- nrow(boot.t)
    neval <- ncol(boot.t)
    if (band.type == "pointwise") {
      probs <- c(alpha / 2.0, 1.0 - alpha / 2.0)
      return(t(apply(boot.t, 2, quantile, probs = probs)))
    }
    if (band.type == "bonferroni") {
      probs <- c(alpha / (2.0 * neval), 1.0 - alpha / (2.0 * neval))
      return(t(apply(boot.t, 2, quantile, probs = probs)))
    }
    if (band.type == "simultaneous") {
      return(old_scs(boot.t, conf.level = 1.0 - alpha)$conf.int)
    }
    if (band.type == "all") {
      return(list(
        pointwise = old_bounds(boot.t, alpha, "pointwise"),
        bonferroni = old_bounds(boot.t, alpha, "bonferroni"),
        simultaneous = old_bounds(boot.t, alpha, "simultaneous")
      ))
    }
    stop("unsupported band type")
  }

  set.seed(20260310L)
  boot.t <- matrix(sample(c(-2, -1, 0, 0, 1, 2, 3), size = 257L * 9L, replace = TRUE), nrow = 257L, ncol = 9L)
  alpha <- 0.1

  expect_equal(new_bounds(boot.t, alpha, "pointwise", warn.coverage = FALSE), old_bounds(boot.t, alpha, "pointwise"))
  expect_equal(new_bounds(boot.t, alpha, "bonferroni", warn.coverage = FALSE), old_bounds(boot.t, alpha, "bonferroni"))
  expect_equal(new_bounds(boot.t, alpha, "simultaneous", warn.coverage = FALSE), old_bounds(boot.t, alpha, "simultaneous"))

  all.new <- new_bounds(boot.t, alpha, "all", warn.coverage = FALSE)
  all.old <- old_bounds(boot.t, alpha, "all")
  expect_equal(all.new$pointwise, all.old$pointwise)
  expect_equal(all.new$bonferroni, all.old$bonferroni)
  expect_equal(all.new$simultaneous, all.old$simultaneous)

  expect_equal(new_scs(boot.t, conf.level = 0.9, alternative = "two.sided")$conf.int,
               old_scs(boot.t, conf.level = 0.9, alternative = "two.sided")$conf.int)
  expect_equal(new_scs(boot.t, conf.level = 0.9, alternative = "less")$conf.int,
               old_scs(boot.t, conf.level = 0.9, alternative = "less")$conf.int)
  expect_equal(new_scs(boot.t, conf.level = 0.9, alternative = "greater")$conf.int,
               old_scs(boot.t, conf.level = 0.9, alternative = "greater")$conf.int)
})
