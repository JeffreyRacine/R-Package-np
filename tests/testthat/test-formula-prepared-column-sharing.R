test_that("prepared formula columns preserve exact model-frame values", {
  frame <- getFromNamespace(".np_formula_model_frame", "np")
  i <- seq_len(24L)
  d <- data.frame(y = cos(i / 3), x = sin(i / 2),
                  g = factor(rep(letters[1:3], 8)))
  calls <- 0L
  policy <- function(x) { calls <<- calls + 1L; stats::na.omit(x) }
  original <- d
  expected <- stats::model.frame(y ~ x + g, data = d, na.action = stats::na.omit)
  actual <- frame(y ~ x + g, data = d, na.action = policy)
  expect_identical(calls, 1L)
  expect_identical(names(attributes(actual)), names(attributes(expected)))
  expect_identical(d, original)
  expect_identical(actual$x, d$x)
  expect_identical(actual$g, d$g)
  expect_identical(.row_names_info(actual, 0L), .row_names_info(expected, 0L))
  if (capabilities("profmem")) {
    expect_identical(tracemem(actual$x), tracemem(d$x))
    expect_identical(tracemem(actual$g), tracemem(d$g))
    untracemem(actual$x)
    untracemem(actual$g)
  }
  actual$x[1L] <- 999
  expect_identical(d, original)
  d$g[2L] <- "a"
  expect_identical(actual$g, original$g)
})

test_that("prepared sharing retains policy changes and row ownership", {
  frame <- getFromNamespace(".np_formula_model_frame", "np")
  i <- seq_len(24L)
  d <- data.frame(y = cos(i / 3), x = sin(i / 2),
                  g = factor(rep(letters[1:3], 8)))
  for (kind in c("modified", "reordered", "omitted", "custom")) {
    calls <- 0L
    policy <- function(x) {
      calls <<- calls + 1L
      if (kind == "modified") x$x <- x$x + 1
      if (kind == "reordered") x <- x[24:1, , drop = FALSE]
      if (kind == "omitted") x <- x[-1L, , drop = FALSE]
      if (kind == "custom") attr(x$x, "source") <- "NA policy"
      x
    }
    expected <- stats::model.frame(y ~ x + g, data = d, na.action = policy)
    calls <- 0L
    actual <- frame(y ~ x + g, data = d, na.action = policy)
    expect_identical(calls, 1L)
    expect_identical(names(attributes(actual)), names(attributes(expected)))
    expect_identical(actual$x, expected$x)
    expect_identical(actual$g, expected$g)
    expect_identical(.row_names_info(actual, 0L),
                     .row_names_info(expected, 0L))
  }
})

test_that("prepared sharing adds no frame or factor method dispatch", {
  share <- getFromNamespace(".np_formula_share_prepared_columns", "np")
  x <- data.frame(y = cos(1:12), g = factor(rep(letters[1:3], 4)))
  a <- list(x$y, x$g)
  registry <- get(".__S3MethodsTable__.", baseenv())
  key <- "length.data.frame"
  old <- get0(key, envir = registry, inherits = FALSE)
  on.exit(if (is.null(old)) rm(list = key, envir = registry) else
            assign(key, old, envir = registry), add = TRUE)
  hits <- 0L
  registerS3method("length", "data.frame", function(x) {
    hits <<- hits + 1L
    length(unclass(x))
  }, envir = baseenv())
  actual <- share(x, a)
  expect_identical(hits, 0L)
  expect_identical(actual, x)
  expect_identical(names(attributes(actual)), names(attributes(x)))
})


test_that("prepared sharing leaves S4-marked owners intact", {
  share <- getFromNamespace(".np_formula_share_prepared_columns", "np")
  x <- asS4(data.frame(x = sin(1:12), g = factor(rep(letters[1:3], 4))))
  expect_identical(share(x, list(x$x, x$g)), x)
  x <- data.frame(x = sin(1:12), g = factor(rep(letters[1:3], 4)))
  a <- asS4(list(x$x, x$g))
  expect_identical(share(x, a), x)
})

test_that("prepared sharing preserves factor-level encoding marks", {
  share <- getFromNamespace(".np_formula_share_prepared_columns", "np")
  utf8 <- intToUtf8(233L)
  latin1 <- iconv(utf8, from = "UTF-8", to = "latin1")
  x <- data.frame(g = factor(rep(latin1, 4L)))
  aligned <- list(factor(rep(utf8, 4L)))
  actual <- share(x, aligned)
  expect_identical(Encoding(levels(actual$g)), Encoding(levels(x$g)))
  expect_false(identical(Encoding(levels(actual$g)),
                         Encoding(levels(aligned[[1L]]))))
})

test_that("prepared sharing preserves compact sequence representation", {
  frame <- getFromNamespace(".np_formula_model_frame", "np")
  d <- data.frame(k = seq_len(2048L), g = factor(rep(letters[1:2], 1024L)))
  actual <- frame(~ k + g, data = d, na.action = stats::na.omit)
  original <- serialize(d$k, NULL)
  copied <- serialize(d$k[rep(TRUE, nrow(d))], NULL)
  expect_identical(actual$k, d$k)
  expect_identical(serialize(actual$k, NULL), original)
  # Some R builds may choose the same representation for both operations.
  # Where the base copy expands the sequence, it is the negative control.
  if (length(copied) > length(original))
    expect_lt(length(serialize(actual$k, NULL)), length(copied))
})
