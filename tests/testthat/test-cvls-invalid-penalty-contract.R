.cvls_penalty_next_up <- function(value) {
  min_subnormal <- .Machine$double.xmin * .Machine$double.eps
  if (value == 0)
    return(min_subnormal)
  bytes <- writeBin(as.double(value), raw(), size = 8L,
                    endian = .Platform$endian)
  if (.Platform$endian != "little")
    bytes <- rev(bytes)
  carry <- 1L
  for (i in seq_along(bytes)) {
    if (!carry)
      break
    next_byte <- as.integer(bytes[[i]]) + carry
    bytes[[i]] <- as.raw(next_byte %% 256L)
    carry <- next_byte >= 256L
  }
  if (.Platform$endian != "little")
    bytes <- rev(bytes)
  readBin(bytes, what = "double", n = 1L, size = 8L,
          endian = .Platform$endian)
}

.cvls_penalty_c_mse <- function(response) {
  total <- 0
  for (value in response)
    total <- total + value
  center <- total / length(response)
  result <- 0
  for (value in response) {
    difference <- value - center
    result <- result + difference * difference
  }
  result / length(response)
}

.cvls_penalty_reference <- function(response) {
  n <- length(response)
  leave_one_scale <- n / (n - 1)
  leave_one_scale * leave_one_scale * .cvls_penalty_c_mse(response)
}

.cvls_penalty_expected <- function(response, multiplier) {
  reference <- .cvls_penalty_reference(response)
  candidate <- multiplier * reference
  if (candidate == reference)
    candidate <- .cvls_penalty_next_up(reference)
  if (!is.finite(candidate) || candidate <= reference ||
      candidate >= .Machine$double.xmax)
    .Machine$double.xmax
  else
    candidate
}

.cvls_penalty_bad_bandwidth <- function(x, y) {
  bws <- npregbw(
    xdat = x, ydat = y, bws = 1, bandwidth.compute = FALSE,
    bwtype = "fixed", bwscaling = FALSE, ckertype = "uniform",
    bwmethod = "cv.ls", regtype = "lc"
  )
  bws$bw[] <- 0
  bws
}

test_that("CVLS invalid guidance uses the exact delete-one intercept null", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  fixtures <- list(
    n2 = c(0, 1),
    ordinary = c(-2, -1, 0, 0.5, 1, 2, 4, 8),
    constant = rep(7, 8),
    translated = c(-2, -1, 0, 0.5, 1, 2, 4, 8) + 1e6
  )

  for (fixture in names(fixtures)) {
    y <- fixtures[[fixture]]
    x <- data.frame(x = seq(-1, 1, length.out = length(y)))
    bad <- .cvls_penalty_bad_bandwidth(x, y)
    for (multiplier in c(1, 10)) {
      value <- np:::.npregbw_eval_only(
        xdat = x, ydat = y, bws = bad,
        invalid.penalty = "baseline",
        penalty.multiplier = multiplier
      )$objective
      expect_identical(
        as.double(value),
        as.double(.cvls_penalty_expected(y, multiplier)),
        info = paste(fixture, multiplier)
      )
    }
    dbmax <- np:::.npregbw_eval_only(
      xdat = x, ydat = y, bws = bad,
      invalid.penalty = "dbmax", penalty.multiplier = 10
    )$objective
    expect_identical(as.double(dbmax), .Machine$double.xmax,
                     info = fixture)
  }
})
