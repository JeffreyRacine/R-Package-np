test_that("formula selectors and frames consume data and NA factories once", {
  old <- options(np.messages = FALSE)
  on.exit(options(old))
  d <- data.frame(x = seq(-1, 1, length.out = 24L), y = sin(seq(-1, 1, length.out = 24L)))
  count <- c(data = 0L, policy = 0L, na = 0L)
  producer <- function() { count['data'] <<- count['data'] + 1L; d }
  policy <- function() {
    count['policy'] <<- count['policy'] + 1L
    function(z) { count['na'] <<- count['na'] + 1L; stats::na.omit(z) }
  }
  a <- npreg(formula = y ~ x, data = producer(), bws = .4,
             bandwidth.compute = FALSE, na.action = policy())
  expect_identical(unname(count), rep(1L, 3L))
  b <- npreg(formula = y ~ x, data = d, bws = .4, bandwidth.compute = FALSE)
  expect_identical(fitted(a), fitted(b))
  expect_false(grepl('evalq|<environment', paste(deparse(a$bws$call), collapse = ' ')))
  count[] <- 0L
  bw <- npcdensbw(y ~ x, data = producer(), bws = c(.4, .4),
                 bandwidth.compute = FALSE, na.action = policy())
  expect_identical(unname(count), rep(1L, 3L))
  expect_identical(bw$call$data, quote(producer()))
  before <- count
  fit <- npcdens(bw)
  copy <- unserialize(serialize(bw, NULL))
  expect_equal(fitted(fit), fitted(npcdens(copy)), tolerance = 0)
  expect_identical(count, before)
  count[] <- 0L
  sums <- npksum(~x, data = producer(), bws = .4, na.action = policy())
  expect_identical(unname(count), rep(1L, 3L))
})

test_that("CDF formula grid callbacks are not replayed after grid admission", {
  old <- options(np.messages = FALSE)
  on.exit(options(old))
  d <- data.frame(x = seq(-1, 1, length.out = 24L), y = seq(0, 1, length.out = 24L))
  count <- 0L
  grid <- function() { count <<- count + 1L; d[1:5, , drop = FALSE] }
  a <- npudistbw(~x, data = d, gdata = grid(), bws = .4, bandwidth.compute = FALSE)
  expect_identical(count, 1L)
  expect_identical(a$call$gdata, quote(grid()))
  count <- 0L
  b <- npcdistbw(y ~ x, data = d, gdata = grid(), bws = c(.4, .4), bandwidth.compute = FALSE)
  expect_identical(count, 1L)
  expect_identical(b$call$gdata, quote(grid()))
})
