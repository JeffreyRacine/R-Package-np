test_that("prepared formula values retain ordinary column sharing", {
  skip_if_not(capabilities("profmem"))
  d <- data.frame(x = sin(seq_len(128L)), y = cos(seq_len(128L)))
  before <- vapply(d, tracemem, character(1L))
  on.exit(lapply(d, untracemem), add = TRUE)
  frame <- .np_formula_model_frame(y ~ x, data = d, na.action = stats::na.omit)
  after <- vapply(frame, tracemem, character(1L))
  on.exit(lapply(frame, untracemem), add = TRUE)
  expect_identical(after, before[names(frame)])
  old <- frame$x
  d$x[1L] <- 100
  expect_identical(frame$x, old)
  expect_false(identical(frame$x, d$x))
})

test_that("prepared values preserve callback and row semantics under quote masking", {
  d <- data.frame(x = sin(seq_len(16L)), y = cos(seq_len(16L)))
  d$x[3L] <- NA_real_
  hits <- c(transform = 0L, action = 0L)
  shift <- function(x) { hits[["transform"]] <<- hits[["transform"]] + 1L; x + 2 }
  policy <- function(x) { hits[["action"]] <<- hits[["action"]] + 1L; stats::na.exclude(x) }
  quote <- function(...) stop("caller quote must not own prepared values")
  f <- y ~ shift(x)
  frame <- .np_formula_model_frame(f, data = d, subset = seq_len(12L), na.action = policy)
  keep <- setdiff(seq_len(12L), 3L)
  expect_identical(unname(hits), c(1L, 1L))
  expect_identical(frame[[1L]], d$y[keep])
  expect_identical(frame[[2L]], (d$x + 2)[keep])
  expect_identical(row.names(frame), as.character(keep))
  expect_s3_class(attr(frame, "na.action"), "exclude")
  expect_identical(as.integer(attr(frame, "na.action")), 3L)
})
