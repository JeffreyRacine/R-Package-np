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

test_that("stored formula reentry honors explicit data override", {
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

  bw.cm <- npcdensbw(
    yc ~ x,
    data = data_a,
    bws = c(0.45, 0.45),
    bandwidth.compute = FALSE
  )
  .formula_reentry_expect(
    npconmode(bws = bw.cm, data = data_b, probabilities = TRUE)$probabilities,
    npconmode(bws = bw.cm, txdat = data_b["x"], tydat = data_b$yc,
              probabilities = TRUE)$probabilities,
    npconmode(bws = bw.cm, probabilities = TRUE)$probabilities,
    npconmode(bws = bw.cm, data = NULL, probabilities = TRUE)$probabilities
  )
  .formula_reentry_expect(
    npconmode(bws = bw.cm, data = data_b, newdata = eval_b["x"],
              probabilities = TRUE)$probabilities,
    npconmode(bws = bw.cm, txdat = data_b["x"], tydat = data_b$yc,
              exdat = eval_b["x"], probabilities = TRUE)$probabilities,
    npconmode(bws = bw.cm, newdata = eval_b["x"],
              probabilities = TRUE)$probabilities,
    npconmode(bws = bw.cm, data = NULL, newdata = eval_b["x"],
              probabilities = TRUE)$probabilities
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

  bw.sc <- npscoefbw(
    y ~ x | z,
    data = data_a,
    bws = 0.45,
    bandwidth.compute = FALSE
  )
  .formula_reentry_expect(
    npscoef(bws = bw.sc, data = data_b, errors = FALSE)$mean,
    npscoef(bws = bw.sc, txdat = data_b["x"], tzdat = data_b["z"],
            tydat = data_b$y, errors = FALSE)$mean,
    npscoef(bws = bw.sc, errors = FALSE)$mean,
    npscoef(bws = bw.sc, data = NULL, errors = FALSE)$mean
  )
  .formula_reentry_expect(
    npscoef(bws = bw.sc, data = data_b, newdata = eval_b,
            errors = FALSE)$mean,
    npscoef(bws = bw.sc, txdat = data_b["x"], tzdat = data_b["z"],
            tydat = data_b$y, exdat = eval_b["x"], ezdat = eval_b["z"],
            errors = FALSE)$mean,
    npscoef(bws = bw.sc, newdata = eval_b, errors = FALSE)$mean,
    npscoef(bws = bw.sc, data = NULL, newdata = eval_b,
            errors = FALSE)$mean
  )

  bw.pl <- npplregbw(
    y ~ z | x,
    data = data_a,
    bws = matrix(c(0.45, 0.45), nrow = 2L),
    bandwidth.compute = FALSE
  )
  .formula_reentry_expect(
    npplreg(bws = bw.pl, data = data_b)$mean,
    npplreg(bws = bw.pl, txdat = data_b["z"], tzdat = data_b["x"],
            tydat = data_b$y)$mean,
    npplreg(bws = bw.pl)$mean,
    npplreg(bws = bw.pl, data = NULL)$mean
  )
  .formula_reentry_expect(
    npplreg(bws = bw.pl, data = data_b, newdata = eval_b)$mean,
    npplreg(bws = bw.pl, txdat = data_b["z"], tzdat = data_b["x"],
            tydat = data_b$y, exdat = eval_b["z"],
            ezdat = eval_b["x"])$mean,
    npplreg(bws = bw.pl, newdata = eval_b)$mean,
    npplreg(bws = bw.pl, data = NULL, newdata = eval_b)$mean
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
