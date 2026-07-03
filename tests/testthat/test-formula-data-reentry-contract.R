.formula_reentry_data <- function(n = 36L) {
  x <- seq(-1, 1, length.out = n)
  z <- seq(0.15, 1.15, length.out = n)
  data_a <- data.frame(
    x = x,
    z = z,
    y = sin(2 * x) + 0.35 * z,
    yc = factor(ifelse(x + 0.1 * z > 0, "b", "a"), levels = c("a", "b"))
  )
  data_b <- data.frame(
    x = rev(x) + 0.17,
    z = rev(z) - 0.09,
    y = -0.6 + cos(2 * rev(x)) - 0.45 * rev(z),
    yc = factor(ifelse(rev(x) - 0.15 * rev(z) > -0.1, "b", "a"),
                levels = c("a", "b"))
  )
  eval_b <- data.frame(
    x = seq(-0.65, 0.75, length.out = 7L),
    z = seq(0.25, 0.95, length.out = 7L)
  )
  list(data_a = data_a, data_b = data_b, eval_b = eval_b)
}

.formula_reentry_flat <- function(x) {
  unname(as.numeric(x))
}

.formula_reentry_expect <- function(reentry, native, stale, null,
                                    tolerance = 1e-10) {
  reentry <- .formula_reentry_flat(reentry)
  native <- .formula_reentry_flat(native)
  stale <- .formula_reentry_flat(stale)
  null <- .formula_reentry_flat(null)

  expect_equal(reentry, native, tolerance = tolerance)
  expect_equal(null, stale, tolerance = tolerance)
  expect_gt(max(abs(reentry - stale), na.rm = TRUE), 1e-7)
}

.formula_reentry_capture_default <- function(method, expr) {
  expr <- substitute(expr)
  capture <- new.env(parent = emptyenv())
  assign(".formula_reentry_capture", capture, envir = .GlobalEnv)
  on.exit(rm(".formula_reentry_capture", envir = .GlobalEnv), add = TRUE)

  trace(
    method,
    tracer = quote({
      cap <- get(".formula_reentry_capture", envir = .GlobalEnv)
      if (exists("txdat", inherits = FALSE))
        cap$txdat <- txdat
      if (exists("tydat", inherits = FALSE))
        cap$tydat <- tydat
      if (exists("tzdat", inherits = FALSE))
        cap$tzdat <- tzdat
      stop("FORMULA_REENTRY_CAPTURE", call. = FALSE)
    }),
    where = asNamespace("npRmpi"),
    print = FALSE
  )
  on.exit(untrace(method, where = asNamespace("npRmpi")), add = TRUE)

  expect_error(eval(expr, parent.frame()), "FORMULA_REENTRY_CAPTURE",
               fixed = TRUE)
  as.list(capture)
}

test_that("stored formula reentry honors explicit data override in executable routes", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old.opts <- options(npRmpi.autodispatch = FALSE, np.messages = FALSE)
  on.exit(options(old.opts), add = TRUE)

  d <- .formula_reentry_data()
  data_a <- d$data_a
  data_b <- d$data_b
  eval_b <- d$eval_b

  bw.qr <- npcdistbw(
    y ~ x,
    data = data_a,
    bws = c(0.45, 0.45),
    bandwidth.compute = FALSE
  )
  .formula_reentry_expect(
    fitted(npqreg(bws = bw.qr, data = data_b, tau = 0.5)),
    fitted(npqreg(bws = bw.qr, txdat = data_b["x"], tydat = data_b$y,
                  tau = 0.5)),
    fitted(npqreg(bws = bw.qr, tau = 0.5)),
    fitted(npqreg(bws = bw.qr, data = NULL, tau = 0.5))
  )
  .formula_reentry_expect(
    fitted(npqreg(bws = bw.qr, data = data_b, newdata = eval_b["x"],
                  tau = 0.5)),
    fitted(npqreg(bws = bw.qr, txdat = data_b["x"], tydat = data_b$y,
                  exdat = eval_b["x"], tau = 0.5)),
    fitted(npqreg(bws = bw.qr, newdata = eval_b["x"], tau = 0.5)),
    fitted(npqreg(bws = bw.qr, data = NULL, newdata = eval_b["x"],
                  tau = 0.5))
  )

  bw.reg <- npregbw(
    y ~ x,
    data = data_a,
    bws = 0.45,
    bandwidth.compute = FALSE
  )
  .formula_reentry_expect(
    npreghat(bws = bw.reg, data = data_b, output = "apply"),
    npreghat(bws = bw.reg, txdat = data_b["x"], y = data_b$y,
             output = "apply"),
    npreghat(bws = bw.reg, output = "apply"),
    npreghat(bws = bw.reg, data = NULL, output = "apply")
  )
  .formula_reentry_expect(
    npreghat(bws = bw.reg, data = data_b, newdata = eval_b["x"],
             output = "apply"),
    npreghat(bws = bw.reg, txdat = data_b["x"], y = data_b$y,
             exdat = eval_b["x"], output = "apply"),
    npreghat(bws = bw.reg, newdata = eval_b["x"], output = "apply"),
    npreghat(bws = bw.reg, data = NULL, newdata = eval_b["x"],
             output = "apply")
  )

  bw.sig <- npregbw(
    y ~ x + z,
    data = data_a,
    bws = c(0.45, 0.45),
    bandwidth.compute = FALSE
  )
  .formula_reentry_expect(
    npsigtest(bws = bw.sig, data = data_b, boot.num = 9L,
              random.seed = 20260703L)$In,
    npsigtest(bws = bw.sig, xdat = data_b[c("x", "z")], ydat = data_b$y,
              boot.num = 9L, random.seed = 20260703L)$In,
    npsigtest(bws = bw.sig, boot.num = 9L, random.seed = 20260703L)$In,
    npsigtest(bws = bw.sig, data = NULL, boot.num = 9L,
              random.seed = 20260703L)$In
  )
})

test_that("stored formula reentry passes explicit data to stalled MPI fit routes", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old.opts <- options(npRmpi.autodispatch = FALSE, np.messages = FALSE)
  on.exit(options(old.opts), add = TRUE)

  d <- .formula_reentry_data()
  data_a <- d$data_a
  data_b <- d$data_b

  bw.cm <- npcdensbw(
    yc ~ x,
    data = data_a,
    bws = c(0.45, 0.45),
    bandwidth.compute = FALSE
  )
  cap.cm <- .formula_reentry_capture_default(
    "npconmode.conbandwidth",
    npconmode(bws = bw.cm, data = data_b)
  )
  expect_equal(as.numeric(cap.cm$txdat$x), data_b$x, tolerance = 0)
  expect_identical(as.character(cap.cm$tydat$yc), as.character(data_b$yc))

  bw.sc <- npscoefbw(
    y ~ x | z,
    data = data_a,
    bws = 0.45,
    bandwidth.compute = FALSE
  )
  cap.sc <- .formula_reentry_capture_default(
    "npscoef.scbandwidth",
    npscoef(bws = bw.sc, data = data_b, errors = FALSE)
  )
  expect_equal(as.numeric(cap.sc$txdat$x), data_b$x, tolerance = 0)
  expect_equal(as.numeric(cap.sc$tzdat$z), data_b$z, tolerance = 0)
  expect_equal(as.numeric(cap.sc$tydat), data_b$y, tolerance = 0)

  bw.pl <- npplregbw(
    y ~ z | x,
    data = data_a,
    bws = matrix(c(0.45, 0.45), nrow = 2L),
    bandwidth.compute = FALSE
  )
  cap.pl <- .formula_reentry_capture_default(
    "npplreg.plbandwidth",
    npplreg(bws = bw.pl, data = data_b)
  )
  expect_equal(as.numeric(cap.pl$txdat$z), data_b$z, tolerance = 0)
  expect_equal(as.numeric(cap.pl$tzdat$x), data_b$x, tolerance = 0)
  expect_equal(as.numeric(cap.pl$tydat), data_b$y, tolerance = 0)
})
