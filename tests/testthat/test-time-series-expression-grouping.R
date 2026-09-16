test_that("time alignment preserves expression grouping and matrix shape", {
  align <- getFromNamespace(".np_formula_align_values", "np")
  y <- ts(seq_len(12), start = c(2000, 1), frequency = 4)
  x <- ts(cbind(a = 101:112, b = 201:212), start = c(2000, 1), frequency = 4)
  z <- ts(cbind(c = 301:312, d = 401:412, e = 501:512),
          start = c(2000, 1), frequency = 4)
  out <- NULL
  expect_warning(out <- align(list(y = y, x = lag(x, -1), yl = lag(y, -1))), NA)
  expect_identical(names(out), c("y", "x", "yl"))
  expect_identical(vapply(out, NCOL, integer(1)), c(y = 1L, x = 2L, yl = 1L))
  expect_equal(out$y, as.numeric(y)[2:12])
  expect_equal(out$x, unclass(x)[1:11, , drop = FALSE], ignore_attr = TRUE)
  expect_identical(colnames(out$x), c("a", "b"))
  expect_equal(out$yl, as.numeric(y)[1:11])
  expect_false(any(vapply(out, inherits, logical(1), "ts")))
  two <- align(list(x = x, z = lag(z, -2)))
  expect_equal(two$x, unclass(x)[3:12, , drop = FALSE], ignore_attr = TRUE)
  expect_equal(two$z, unclass(z)[1:10, , drop = FALSE], ignore_attr = TRUE)
  expect_identical(colnames(two$z), c("c", "d", "e"))
  one <- ts(matrix(21:32, ncol = 1, dimnames = list(NULL, "only")),
            start = c(2000, 1), frequency = 4)
  singleton <- align(list(y = y, one = lag(one, -1)))
  expect_identical(dim(singleton$one), c(11L, 1L))
  expect_identical(colnames(singleton$one), "only")
  expect_equal(as.numeric(singleton$one), as.numeric(one)[1:11])
  lone <- align(list(x = x))
  expect_identical(dim(lone$x), c(12L, 2L))
  expect_identical(colnames(lone$x), c("a", "b"))
})

test_that("time alignment leaves ordinary values and subset inputs intact", {
  align <- getFromNamespace(".np_formula_align_values", "np")
  ordinary <- list(y = c(1, NA, 3), f = ordered(c("b", "a", "b")),
                   m = matrix(1:6, ncol = 2))
  expect_identical(align(ordinary), ordinary)
  y <- ts(1:12, start = c(2001, 1), frequency = 12)
  x <- ts(cbind(a = 21:32, b = 41:52), start = c(2001, 1), frequency = 12)
  x[3, 2] <- NA
  f <- factor(rep(c("a", "b"), 5))
  out <- align(list(f = f, x = x, y = lag(y, -2)))
  expect_identical(out$f, f)
  expect_identical(dim(out$x), c(10L, 2L))
  expect_true(is.na(out$x[1, 2]))
  expect_equal(as.numeric(out$y), 1:10)
  expect_error(suppressWarnings(align(list(y = y,
    z = ts(1:3, start = 2010, frequency = 12)))), "no common observations")
})

test_that("aligned frames retain matrix variables without inventing family support", {
  align.terms <- getFromNamespace(".np_formula_aligned_terms", "np")
  y <- ts(seq_len(12), start = c(2000, 1), frequency = 4)
  x <- ts(cbind(a = 101:112, b = 201:212), start = c(2000, 1), frequency = 4)
  f <- y ~ lag(x, -1) + lag(y, -1)
  mf <- model.frame(align.terms(terms(f)))
  expect_identical(nrow(mf), 11L)
  expect_equal(unname(model.response(mf)), as.numeric(y)[2:12])
  expect_identical(dim(mf[[2L]]), c(11L, 2L))
  expect_equal(mf[[3L]], as.numeric(y)[1:11])
  # These families did not support a matrix stored as one formula variable.
  # In particular, never accept its first column as if it were the variable.
  expect_warning(expect_error(npindex(y ~ lag(x, -1), bws = c(1, .8), se = FALSE),
    "xdat must contain at least one continuous variable"), NA)
  expect_warning(expect_error(npksum(f, bws = c(.8, .8)),
    "supplied bandwidths do not match 'txdat' in type"), NA)
})
