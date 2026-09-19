test_that("scalar numeric AsIs formula terms match ordinary numeric transforms", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(191013)
  d <- data.frame(x = runif(30, -1, 1), y = rnorm(30))
  square <- function(x) x^2
  for (family in c("npreg", "npindex", "npcdens", "npcdist", "npudens", "npudist")) {
    h <- switch(family, npindex = c(1, .6), npcdens = c(.7, .6),
                npcdist = c(.7, .6), .6)
    f <- if (family %in% c("npudens", "npudist")) ~I(x^2) else y ~ I(x^2)
    g <- if (family %in% c("npudens", "npudist")) ~square(x) else y ~ square(x)
    constructor <- get(paste0(family, "bw"))
    a <- constructor(formula = f, data = d, bws = h, bandwidth.compute = FALSE)
    b <- constructor(formula = g, data = d, bws = h, bandwidth.compute = FALSE)
    af <- get(family)(a, se = FALSE)
    bf <- get(family)(b, se = FALSE)
    expect_equal(fitted(af), fitted(bf), tolerance = 1e-12)
    expect_equal(predict(af, newdata = d[1:4, ]),
                 predict(bf, newdata = d[1:4, ]), tolerance = 1e-12)
    expect_identical(a$formula, f)
  }
  count <- 0L
  once <- function(x) { count <<- count + 1L; I(x^2) }
  npreg(y ~ once(x), data = d, bws = .6, bandwidth.compute = FALSE)
  expect_identical(count, 1L)
  expect_error(npreg(y ~ I(cbind(x, x)), data = d, bws = .6,
                     bandwidth.compute = FALSE), "matrix|in type")
})
