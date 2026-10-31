test_that("public regression checks selected outcomes without changing predictors", {
  withr::local_options(np.messages = FALSE)
  n <- 24L
  x <- data.frame(x = seq(.05, .95, length.out = n),
                  f = factor(rep(c("a", "b", "c"), 8)),
                  o = ordered(rep(c(10, 30, 50), 8)))
  y <- sin(x$x)
  b <- npregbw(xdat = x, ydat = y, bws = c(.3, .2, .2),
               bandwidth.compute = FALSE)
  bad <- list(x$f, x$o, rep(TRUE, n), as.character(y), as.complex(y),
              as.raw(seq_len(n)), as.list(y), NULL, matrix(y, n))
  for (response in bad) {
    expect_error(npregbw(xdat = x, ydat = response, bws = b,
                         bandwidth.compute = FALSE), "must be a numeric")
    expect_error(npreg(b, txdat = x, tydat = response), "must be a numeric")
    expect_error(npreg(b, txdat = x, tydat = y, eydat = response),
                 "must be a numeric")
    expect_length(fitted(npreg(b)), n)
  }
  # A private context cannot be supplied through public dots to grant an exemption.
  expect_error(npreg(b, txdat = x, tydat = x$f,
                      .np.context = list(public = FALSE)))
  expect_error(npregbw(xdat = x, ydat = x$f, bws = c(.3,.2,.2),
                        bandwidth.compute = FALSE,
                        .np.context = list(public = FALSE)))
  calls <- 0L
  response <- function() { calls <<- calls + 1L; y }
  numeric <- npregbw(xdat = x, ydat = response(), bws = c(.3,.2,.2),
                     bandwidth.compute = FALSE)
  expect_identical(calls, 1L)
  expect_length(fitted(npreg(numeric)), n)
  expect_identical(calls, 1L)
})

test_that("categorical children retain their scales only inside auxiliary and operator roles", {
  withr::local_options(np.messages = FALSE)
  n <- 24L
  z <- data.frame(z = seq(.05, .95, length.out = n))
  y <- sin(seq_len(n)/3)
  cats <- list(factor(rep(c("30", "10", "50"), 8),
                      levels = c("50", "10", "30", "70")),
               ordered(rep(c("30", "10", "50"), 8),
                       levels = c("10", "30", "50", "70")))
  scores <- list(rep(c(3,2,1), 8), rep(c(30,10,50), 8))
  for (j in seq_along(cats)) {
    x <- data.frame(cat = cats[[j]])
    b <- npplregbw(xdat = x, ydat = y, zdat = z,
                   bws = matrix(.3, 2, 1), bandwidth.compute = FALSE)
    reference <- npplregbw(xdat = data.frame(cat = scores[[j]]), ydat = y,
                   zdat = z, bws = matrix(.3, 2, 1), bandwidth.compute = FALSE)
    expect_identical(fitted(npplreg(b)), fitted(npplreg(reference)))
    # Informative parent progress keeps child fits at the existing public MPI
    # entry; it must use the same scores and retain the original child object.
    progress.fit <- withr::with_options(list(np.messages = TRUE), npplreg(b))
    expect_identical(fitted(progress.fit), fitted(npplreg(reference)))
    expect_identical(progress.fit$bws$bw[[2L]], b$bw[[2L]])
    child <- unserialize(serialize(b$bw[[2L]], NULL))
    expect_identical(child[[".np.native.training"]]$ydat, cats[[j]])
    expect_error(npreg(child), "must be a numeric")
    expect_error(npregbw(xdat = z, ydat = cats[[j]], bws = child,
                         bandwidth.compute = FALSE), "must be a numeric")
    H <- npreghat(child, txdat = z)
    expect_equal(as.numeric(npreghat(child, output = "apply")),
                  as.numeric(H %*% scores[[j]]), tolerance = 2e-11)
  }
})


test_that("regression forwarding preserves public errors and signals them once", {
  withr::local_options(np.messages = FALSE)
  x <- data.frame(x = seq(0, 1, length.out = 24))
  y <- sin(x$x)
  b <- npregbw(xdat = x, ydat = y, bws = .3, bandwidth.compute = FALSE)
  seen <- 0L
  error <- tryCatch(withCallingHandlers(
    npreg(b, txdat = x, tydat = y[-1L]),
    error = function(condition) seen <<- seen + 1L), error = identity)
  expect_identical(seen, 1L)
  expect_match(conditionMessage(error), "number of explanatory data")
  # MPI retains its existing opcode condition owner; local errors retain the
  # public S3 call and never expose the invocation-local forwarding context.
  expect_false(any(grepl(".np.context", deparse(conditionCall(error)), fixed = TRUE)))
  expect_false(is.function(conditionCall(error)[[1L]]))
  expect_length(fitted(npreg(b)), nrow(x))
})
