t1_npsig_fixture <- function() {
  n <- 20L
  xdat <- data.frame(
    z = factor(rep(c("a", "b"), length.out = n)),
    x = seq(-1, 1, length.out = n)
  )
  ydat <- 0.35 * (xdat$z == "b") + sin(1.3 * xdat$x) +
    seq_len(n) / 1000
  bws <- npregbw(
    xdat = xdat,
    ydat = ydat,
    bws = c(0.2, 0.55),
    bandwidth.compute = FALSE,
    regtype = "ll",
    bwtype = "fixed"
  )
  list(xdat = xdat, ydat = ydat, bws = bws)
}

t1_npsig_rademacher_oracle <- function(fixture, B, joint, seed, index = 1:2) {
  pkg <- environmentName(environment(npsigtest))
  streamed.joint <- getFromNamespace(".np_npsig_streamed_response_statistic", pkg)
  streamed.individual <- getFromNamespace(".np_npsig_streamed_iid_tile", pkg)
  cast <- getFromNamespace("cast", pkg)
  xdat <- fixture$xdat
  ydat <- fixture$ydat
  bws <- fixture$bws
  n <- nrow(xdat)

  withr::with_seed(seed, {
    unrestricted.residuals <- scale(residuals(npreg(bws = bws)))
    unrestricted.scale <- attr(unrestricted.residuals, "scaled:scale")
    unrestricted.center <- attr(unrestricted.residuals, "scaled:center")

    make.null <- function(tested) {
      null.frame <- xdat
      for (variable in tested) {
        median.value <- uocquantile(xdat[[variable]], 0.5)
        null.frame[[variable]] <- if (is.factor(xdat[[variable]]) ||
                                      is.ordered(xdat[[variable]])) {
          cast(median.value, xdat[[variable]], same.levels = TRUE)
        } else {
          median.value
        }
      }
      null.mean <- npreg(
        txdat = xdat, tydat = ydat, exdat = null.frame, bws = bws
      )$mean
      residual.pool <- as.numeric(
        scale(ydat - null.mean) * unrestricted.scale + unrestricted.center
      )
      list(
        mean = null.mean,
        residuals = residual.pool - mean(residual.pool)
      )
    }

    make.responses <- function(null) {
      vapply(seq_len(B), function(unused) {
        uniforms <- stats::runif(n)
        multiplier <- ifelse(uniforms <= 0.5, -1, 1)
        null$mean + null$residuals * multiplier
      }, numeric(n))
    }

    if (isTRUE(joint)) {
      null <- make.null(index)
      responses <- make.responses(null)
      values <- unlist(lapply(seq.int(1L, B, by = 8L), function(start) {
        columns <- start:min(start + 7L, B)
        streamed.joint(
          bws = bws,
          xdat = xdat,
          index = index,
          response.matrix = responses[, columns, drop = FALSE],
          pivotal = FALSE
        )
      }), use.names = FALSE)
      return(matrix(values, ncol = 1L))
    }

    output <- matrix(NA_real_, nrow = B, ncol = length(index))
    for (column in seq_along(index)) {
      tested <- index[[column]]
      null <- make.null(tested)
      responses <- make.responses(null)
      output[, column] <- unlist(lapply(
        seq.int(1L, B, by = 8L),
        function(start) {
          columns <- start:min(start + 7L, B)
          streamed.individual(
            bws = bws,
            xdat = xdat,
            tested.index = tested,
            response.matrix = responses[, columns, drop = FALSE],
            null.mean = null$mean,
            residual.pool = null$residuals,
            pivotal = FALSE
          )
        }
      ), use.names = FALSE)
    }
    output
  })
}

test_that("Rademacher probability is distinct from the Mammen probability", {
  uniforms <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
  multipliers <- ifelse(uniforms <= 0.5, -1, 1)

  expect_identical(unique(multipliers), c(-1, 1))
  expect_equal(mean(multipliers == -1), 0.5)
  expect_equal(mean(multipliers == 1), 0.5)
  expect_equal(mean(multipliers), 0)
  expect_equal(mean((multipliers - mean(multipliers))^2), 1)
  expect_equal(mean(multipliers^2), 1)

  implementation <- paste(deparse(getS3method("npsigtest", "rbandwidth")),
                          collapse = "\n")
  pkg <- environmentName(environment(npsigtest))
  uses <- gregexpr("P.rademacher", implementation, fixed = TRUE)[[1L]]
  expect_match(implementation, "P.rademacher <- 0.5", fixed = TRUE)
  expect_length(uses, if (identical(pkg, "npRmpi")) 9L else 5L)
  expect_false(grepl("c(-1, 1, P.a)", implementation, fixed = TRUE))
  expect_false(grepl("num.obs, -1, 1, P.a", implementation, fixed = TRUE))
})

test_that("seeded joint and individual Rademacher pseudo-outcomes match an independent oracle", {
  fixture <- t1_npsig_fixture()

  for (B in c(9L, 17L)) {
    for (joint in c(TRUE, FALSE)) {
      set.seed(701L + B + as.integer(joint))
      before <- .Random.seed
      actual <- npsigtest(
        bws = fixture$bws,
        xdat = fixture$xdat,
        ydat = fixture$ydat,
        B = B,
        boot.method = "wild-rademacher",
        boot.type = "I",
        pivot = FALSE,
        joint = joint,
        index = seq_len(ncol(fixture$xdat)),
        random.seed = 20260902L
      )
      expect_identical(.Random.seed, before)
      expected <- t1_npsig_rademacher_oracle(
        fixture = fixture,
        B = B,
        joint = joint,
        seed = 20260902L
      )
      expect_equal(actual$In.bootstrap, expected, tolerance = 2e-10,
                   info = paste("B", B, "joint", joint))
    }
  }
})
