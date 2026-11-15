test_that("full row selection preserves base values and row-name records", {
  rows <- getFromNamespace(".np_data_rows", "npRmpi")
  plain <- data.frame(x = sin(seq_len(12)), y = seq_len(12))
  inputs <- list(plain, plain[c(3, 7, 11), , drop = FALSE],
    structure(plain, row.names = 101:112),
    structure(plain, row.names = letters[1:12]),
    transform(plain, y = factor(y)), transform(plain, y = ordered(y)),
    structure(plain, names = c("x", "x")), plain[FALSE, , drop = FALSE],
    plain[, FALSE, drop = FALSE])
  for (x in inputs) {
    keep <- rep(TRUE, .row_names_info(x, 2L))
    expected <- x[keep, , drop = FALSE]
    actual <- rows(x, keep)
    expect_identical(actual, expected)
    expect_identical(names(attributes(actual)), names(attributes(expected)))
    expect_identical(.row_names_info(actual, 0L), .row_names_info(expected, 0L))
    expect_identical(as.matrix(actual), as.matrix(expected))
  }
  for (keep in list(rep(c(TRUE, FALSE), 6), 12:1, TRUE,
                    c(NA, rep(TRUE, 11)), setNames(rep(TRUE, 12), letters[1:12])))
    expect_identical(rows(plain, keep), plain[keep, , drop = FALSE])
  for (x in list(seq_len(12), sin(seq_len(12)), setNames(1:12, letters[1:12])))
    expect_identical(rows(x, rep(TRUE, 12), frame = FALSE), x[rep(TRUE, 12)])
})

test_that("row sharing does not change copy-on-modify ownership", {
  rows <- getFromNamespace(".np_data_rows", "npRmpi")
  x <- data.frame(x = sin(seq_len(12)), g = factor(rep(letters[1:3], 4)))
  expected <- x[rep(TRUE, 12), , drop = FALSE]
  retained <- rows(x, rep(TRUE, 12))
  x$x[1L] <- 10
  x$g[2L] <- "a"
  expect_identical(retained, expected)
  retained$x[3L] <- 20
  expect_false(x$x[3L] == 20)
})

test_that("malformed frame records retain base repair behavior", {
  rows <- getFromNamespace(".np_data_rows", "npRmpi")
  x <- data.frame(x = sin(seq_len(12)), y = seq_len(12))
  records <- list(c("a", "a", letters[3:12]), c(NA_character_, letters[2:12]),
                  c(5L, 5L, 3:12), c(NA_integer_, 20L, 3:12),
                  c(NA_integer_, -20L))
  for (record in records) {
    value <- structure(x, row.names = record)
    keep <- rep(TRUE, .row_names_info(value, 2L))
    actual <- rows(value, keep)
    expected <- value[keep, , drop = FALSE]
    expect_identical(actual, expected)
    expect_identical(.row_names_info(actual, 0L), .row_names_info(expected, 0L))
  }
})

test_that("registered data methods retain their extraction effects", {
  rows <- getFromNamespace(".np_data_rows", "npRmpi")
  exercise <- function(generic, class, method, x) {
    registry <- get(".__S3MethodsTable__.", envir = baseenv())
    key <- paste(generic, class, sep = ".")
    old <- get0(key, envir = registry, inherits = FALSE)
    on.exit(if (is.null(old)) rm(list = key, envir = registry) else
              assign(key, old, envir = registry), add = TRUE)
    hits <- 0L
    original <- get("[.data.frame", baseenv())
    environment(method) <- environment()
    registerS3method(generic, class, method, envir = baseenv())
    keep <- rep(TRUE, .row_names_info(x, 2L))
    expected <- x[keep, , drop = FALSE]
    before <- hits
    hits <- 0L
    actual <- rows(x, keep)
    expect_identical(hits, before)
    expect_identical(actual, expected)
  }
  x <- data.frame(x = sin(seq_len(12)), y = seq_len(12))
  exercise("[", "data.frame", function(x, ...) {
    hits <<- hits + 1L; original(x, ...)
  }, x)
  exercise("names", "data.frame", function(x) {
    hits <<- hits + 1L; names(unclass(x))
  }, x)
  exercise("length", "data.frame", function(x) {
    hits <<- hits + 1L; length(unclass(x))
  }, x)
  x$y <- ordered(rep(letters[1:3], 4))
  exercise("[", "ordered", function(x, ...) {
    hits <<- hits + 1L; NextMethod("[")
  }, x)
  exercise("dim", "ordered", function(x) {
    hits <<- hits + 1L; NULL
  }, x)
  x$y <- factor(x$y)
  exercise("dim", "factor", function(x) {
    hits <<- hits + 1L; NULL
  }, x)
})


test_that("workspace methods retain extraction effects", {
  rows <- getFromNamespace(".np_data_rows", "npRmpi")
  exercise <- function(generic, class, method, x) {
    key <- paste(generic, class, sep = ".")
    workspace <- globalenv()
    old <- get0(key, envir = workspace, inherits = FALSE)
    on.exit(if (is.null(old)) rm(list = key, envir = workspace) else
              assign(key, old, envir = workspace), add = TRUE)
    hits <- 0L
    environment(method) <- environment()
    assign(key, method, envir = workspace)
    keep <- rep(TRUE, .row_names_info(x, 2L))
    hits <- 0L
    expected <- x[keep, , drop = FALSE]
    before <- hits
    hits <- 0L
    actual <- rows(x, keep)
    expect_identical(hits, before)
    expect_identical(actual, expected)
    expect_gt(before, 0L)
  }
  x <- data.frame(x = sin(seq_len(12)), y = seq_len(12))
  exercise("names", "data.frame", function(x) {
    hits <<- hits + 1L; names(unclass(x))
  }, x)
  exercise("length", "data.frame", function(x) {
    hits <<- hits + 1L; length(unclass(x))
  }, x)
  x$y <- ordered(rep(letters[1:3], 4))
  exercise("[", "ordered", function(x, ...) {
    hits <<- hits + 1L; NextMethod("[")
  }, x)
  exercise("dim", "ordered", function(x) {
    hits <<- hits + 1L; NULL
  }, x)
  x$y <- factor(x$y)
  exercise("dim", "factor", function(x) {
    hits <<- hits + 1L; NULL
  }, x)
})


test_that("row sharing does not force an active method binding twice", {
  rows <- getFromNamespace(".np_data_rows", "npRmpi")
  exercise <- function(generic, class, method, x) {
    key <- paste0(generic, ".", class)
    workspace <- globalenv()
    present <- exists(key, envir = workspace, inherits = FALSE)
    active <- present && bindingIsActive(key, workspace)
    old <- if (active) activeBindingFunction(key, workspace) else
      if (present) get(key, envir = workspace, inherits = FALSE) else NULL
    if (present) rm(list = key, envir = workspace)
    on.exit({
      rm(list = key, envir = workspace)
      if (active) makeActiveBinding(key, old, workspace) else
        if (present) assign(key, old, envir = workspace)
    }, add = TRUE)
    hits <- 0L
    makeActiveBinding(key, function(value) {
      if (!missing(value)) stop("test method binding is read-only")
      hits <<- hits + 1L
      method
    }, workspace)
    keep <- rep(TRUE, .row_names_info(x, 2L))
    hits <- 0L
    expected <- x[keep, , drop = FALSE]
    before <- hits
    hits <- 0L
    actual <- rows(x, keep)
    expect_gt(before, 0L)
    expect_identical(hits, before)
    expect_identical(actual, expected)
  }
  x <- data.frame(x = sin(seq_len(12)), y = seq_len(12))
  exercise("names", "data.frame", function(x) names(unclass(x)), x)
  exercise("length", "data.frame", function(x) length(unclass(x)), x)
  x$y <- ordered(rep(letters[1:3], 4))
  exercise("[", "ordered", function(x, ...) NextMethod("["), x)
  exercise("dim", "ordered", function(x) NULL, x)
  x$y <- factor(x$y)
  exercise("dim", "factor", function(x) NULL, x)
})


test_that("S4-marked frames and columns retain ordinary extraction", {
  rows <- getFromNamespace(".np_data_rows", "npRmpi")
  plain <- getFromNamespace(".np_data_plain_frame", "npRmpi")
  x <- data.frame(x = sin(seq_len(12)), y = seq_len(12))
  inputs <- list(asS4(x))
  for (g in list(asS4(factor(rep(letters[1:3], 4))),
                 asS4(ordered(rep(letters[1:3], 4))),
                 asS4(sin(seq_len(12))))) {
    d <- list(x = x$x, g = g)
    attr(d, "row.names") <- c(NA_integer_, -12L)
    class(d) <- "data.frame"
    inputs[[length(inputs) + 1L]] <- d
  }
  for (d in inputs) {
    keep <- rep(TRUE, 12L)
    expected <- d[keep, , drop = FALSE]
    actual <- rows(d, keep)
    expect_false(plain(d))
    expect_identical(actual, expected)
    expect_identical(.row_names_info(actual, 0L),
                     .row_names_info(expected, 0L))
  }
})

test_that("S4-marked atomic vectors retain ordinary extraction", {
  rows <- getFromNamespace(".np_data_rows", "npRmpi")
  plain <- getFromNamespace(".np_data_plain_column", "npRmpi")
  x <- asS4(sin(seq_len(12)))
  keep <- rep(TRUE, 12L)
  expect_false(plain(x))
  expect_identical(rows(x, keep, frame = FALSE), x[keep])
})
