manual_ordered_density <- function(dat, bws, okertype) {
  x <- lapply(dat, as.integer)
  support <- lapply(dat, function(z) seq_len(nlevels(z)))
  n <- length(x[[1L]])
  out <- numeric(n)

  ordered_kernel <- function(eval, train, lambda, support, okertype) {
    if (identical(okertype, "liracine")) {
      lambda ^ abs(eval - train) * (1.0 - lambda) / (1.0 + lambda)
    } else if (identical(okertype, "racineliyan")) {
      den <- vapply(train, function(z) sum(lambda ^ abs(z - support)), numeric(1L))
      lambda ^ abs(eval - train) / den
    } else {
      stop("unsupported ordered kernel in test", call. = FALSE)
    }
  }

  for (i in seq_len(n)) {
    w <- rep(1.0, n)
    for (j in seq_along(x)) {
      w <- w * ordered_kernel(x[[j]][i], x[[j]], bws[j], support[[j]], okertype)
    }
    out[i] <- sum(w) / n
  }
  out
}

test_that("npudens categorical profile cache is scoped to each call", {
  options(np.messages = FALSE)

  for (seed in c(20260531L, 20260536L)) {
    set.seed(seed)
    dat <- data.frame(
      o1 = ordered(sample(1:5, 240, TRUE)),
      o2 = ordered(sample(1:4, 240, TRUE))
    )

    for (okertype in c("liracine", "racineliyan")) {
      for (bws in list(c(0.05, 0.10), c(0.25, 0.35), c(0.40, 0.45))) {
        bw <- npudensbw(
          ~ o1 + o2,
          data = dat,
          bws = bws,
          bandwidth.compute = FALSE,
          okertype = okertype
        )
        fit <- npudens(bws = bw)
        expect_equal(
          fitted(fit),
          manual_ordered_density(dat, bws, okertype),
          tolerance = 1e-12,
          info = paste(seed, okertype, paste(bws, collapse = ","))
        )
      }
    }
  }
})
