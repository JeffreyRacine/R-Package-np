r21_cms_fixture <- function() {
  set.seed(4518)
  n <- 43L
  d <- data.frame(x = c(rnorm(33, sd = .2), rnorm(10, mean = 3, sd = 2)))
  d$y <- d$x + 2*sin(d$x) + rt(n, 4)*(1+abs(d$x))
  list(d = d, mean = lm(y ~ x, data = d, x = TRUE, y = TRUE),
       quantile = quantreg::rq(y ~ x, data = d, tau = .5))
}

r21_cms_literal <- function(d, model, quantile, bw, weighted) {
  n <- nrow(d)
  h <- bw$bandwidth$x
  type <- bw$type
  radius <- if(type == "fixed") rep(h, n) else
    vapply(seq_len(n), function(i) sort(abs(d$x[i] - d$x[-i]))[h], numeric(1L))
  K <- matrix(0, n, n)
  for(i in seq_len(n)) for(j in seq_len(n)) if(i != j) {
    r <- radius[if(type == "adaptive_nn") j else i]
    K[i,j] <- dnorm((d$x[i]-d$x[j])/r)/r
  }
  score <- if(quantile) ifelse(residuals(model) <= 0, .5, -.5) else residuals(model)
  fhat <- if(weighted) 1 else rowSums(K)/n
  In <- sum(score * drop(K %*% score)/fhat)/n^2
  Q <- sum(score^2 * drop(K^2 %*% score^2)/fhat^2)/n^2
  list(In = In, Omega.hat = 2*h*Q, Jn = n*In/sqrt(2*Q))
}

test_that("conditional-moment components match literal normalized kernel rows", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  fixture <- r21_cms_fixture()
  d <- fixture$d
  for (family in c("mean", "quantile"))
    for(type in c("fixed", "generalized_nn", "adaptive_nn"))
      for(weighted in c(FALSE, TRUE)) for(scaled in c(FALSE, TRUE)) {
        fun <- if(family == "mean") npcmstest else npqcmstest
        h <- if(type == "fixed") 1.1 else 12
        bw <- npregbw(xdat = d["x"], ydat = d$y, bws = h,
          bandwidth.compute = FALSE, bwtype = type, bwscaling = scaled)
        expected <- r21_cms_literal(d, fixture[[family]], family == "quantile", bw, weighted)
        seed <- .Random.seed
        actual <- do.call(fun, list(xdat = d["x"], ydat = d$y, model = fixture[[family]],
          bws = h, bandwidth.compute = FALSE, bwscaling = scaled,
          bwtype = type, distribution = "asymptotic", density.weighted = weighted))
        for(field in names(expected))
          expect_equal(actual[[field]], expected[[field]], tolerance = 1e-12)
        expect_equal(actual$P, pnorm(expected$Jn, lower.tail = FALSE), tolerance = 1e-12)
        expect_identical(.Random.seed, seed)
      }
})

test_that("conditional-moment bootstrap preserves physical twins and RNG", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  fixture <- r21_cms_fixture()
  d <- fixture$d
  for (family in c("mean", "quantile")) for (type in c("fixed", "generalized_nn"))
    for (pivot in c(FALSE, TRUE)) {
      fun <- if(family == "mean") npcmstest else npqcmstest
      h <- if(type == "fixed") 1.1 else 12
      bw <- npregbw(xdat = d["x"], ydat = d$y, bws = h,
        bandwidth.compute = FALSE, bwtype = type, bwscaling = TRUE)
      args <- list(xdat = d["x"], ydat = d$y, model = fixture[[family]],
        bandwidth.compute = FALSE, bwtype = type, B = 9L,
        density.weighted = FALSE, distribution = "bootstrap", pivot = pivot)
      seed <- .Random.seed
      a <- do.call(fun, c(args, list(bws = h, bwscaling = TRUE)))
      b <- do.call(fun, c(args, list(bws = bw$bandwidth$x, bwscaling = FALSE)))
      for(field in c("In", "Omega.hat", "Jn", "P", "In.bootstrap", "Jn.bootstrap"))
        expect_equal(a[[field]], b[[field]], tolerance = 1e-12)
      expect_identical(.Random.seed, seed)
    }
})
