test_that("npunitest count reduction preserves the summation bootstrap", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)

  x <- c(-1.4, -0.7, -0.1, 0.2, 0.2, 0.8, 1.3, 2.1)
  y <- c(-1.1, -0.5, 0, 0.4, 0.9, 1.2)
  bw.x <- 0.46
  bw.y <- 0.51
  random.seed <- 24681L

  reference.statistic <- function(sample.x, sample.y) {
    density.x <- fitted(npudens(
      tdat = sample.x, edat = sample.x, bws = bw.x
    ))
    density.y <- fitted(npudens(
      tdat = sample.y, edat = sample.x, bws = bw.y
    ))
    0.5 * mean((1 - sqrt(density.y / density.x))^2)
  }

  set.seed(random.seed)
  reference.bootstrap <- numeric(9L)
  for (b in seq_len(9L)) {
    sample.x <- x[sample.int(length(x), length(x), replace = TRUE)]
    sample.y <- x[sample.int(length(x), length(y), replace = TRUE)]
    reference.bootstrap[b] <- reference.statistic(sample.x, sample.y)
  }

  out <- npunitest(
    x,
    y,
    method = "summation",
    bootstrap = TRUE,
    B = 9L,
    bw.x = bw.x,
    bw.y = bw.y,
    random.seed = random.seed
  )

  expect_lte(
    abs(out$Srho - reference.statistic(x, y)),
    64 * .Machine$double.eps
  )
  expect_lte(
    max(abs(out$Srho.bootstrap - reference.bootstrap)),
    64 * .Machine$double.eps
  )
  expect_identical(
    out$P,
    mean(reference.bootstrap > reference.statistic(x, y))
  )
})

test_that("symmetric summation reflection preserves direct densities", {
  helper <- getFromNamespace(
    ".np_entropy_symmetric_gaussian_summation", "np"
  )
  data <- c(-2.1, -0.8, -0.2, -0.2, 0.4, 1.3, 2.6)
  bandwidth <- 0.57
  location <- mean(data)
  rotated <- 2 * location - data
  density <- fitted(npudens(tdat = data, edat = data, bws = bandwidth))
  density.rotated <- fitted(npudens(
    tdat = rotated, edat = data, bws = bandwidth
  ))
  reference <- 0.5 * mean((1 - sqrt(density.rotated / density))^2)

  expect_lte(abs(helper(data, bandwidth) - reference),
             64 * .Machine$double.eps)
})
