library(np)

a1_effect_frames <- function(exdat, index) {
  lower <- upper <- exdat
  value <- exdat[[index]]
  lev <- levels(value)

  if (is.ordered(value)) {
    pos <- match(as.character(value), lev)
    lower[[index]] <- ordered(
      lev[ifelse(pos == 1L, 1L, pos - 1L)],
      levels = lev
    )
    upper[[index]] <- ordered(
      lev[ifelse(pos == 1L, 2L, pos)],
      levels = lev
    )
  } else {
    lower[[index]] <- factor(lev[[1L]], levels = lev)
    upper[[index]] <- factor(as.character(value), levels = lev)
  }

  list(lower = lower, upper = upper)
}

a1_effect_data <- function() {
  set.seed(20260830L)
  n <- 48L
  x <- data.frame(
    z = seq(-1, 1, length.out = n),
    group = factor(
      rep(c("base", "middle", "active"), length.out = n),
      levels = c("base", "middle", "active")
    )
  )
  y <- 0.65 * x$z + 0.35 * (x$group == "middle") +
    0.85 * (x$group == "active") + rnorm(n, sd = 0.11)
  ex <- x[c(2L, 9L, 17L, 26L, 35L, 44L), , drop = FALSE]
  list(x = x, y = y, ex = ex, frames = a1_effect_frames(ex, 2L))
}

test_that("general-LP regression categorical gradients are endpoint effects", {
  dat <- a1_effect_data()
  x.before <- dat$x
  ex.before <- dat$ex
  bw <- npregbw(
    xdat = dat$x,
    ydat = dat$y,
    bws = c(0.42, 0.20),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 2L,
    basis = "tensor"
  )

  fit <- suppressWarnings(npreg(
    bws = bw,
    txdat = dat$x,
    tydat = dat$y,
    exdat = dat$ex,
    gradients = TRUE,
    se = TRUE
  ))
  oracle <- npreg(
    bws = bw,
    txdat = dat$x,
    tydat = dat$y,
    exdat = dat$frames$upper,
    gradients = FALSE
  )$mean - npreg(
    bws = bw,
    txdat = dat$x,
    tydat = dat$y,
    exdat = dat$frames$lower,
    gradients = FALSE
  )$mean

  expect_equal(fit$grad[, 2L], oracle, tolerance = 1e-10)
  expect_true(all(is.finite(fit$gerr[, 2L])))
  expect_true(all(fit$gerr[, 2L] >= 0))
  expect_true(any(fit$gerr[, 2L] > 0))
  expect_identical(dat$x, x.before)
  expect_identical(dat$ex, ex.before)

  ox <- transform(
    dat$x,
    group = ordered(group, levels = levels(group))
  )
  oex <- ox[c(1L, 2L, 3L, 16L, 17L, 18L), , drop = FALSE]
  oframes <- a1_effect_frames(oex, 2L)
  obw <- npregbw(
    xdat = ox,
    ydat = dat$y,
    bws = c(0.42, 0.20),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 2L,
    basis = "tensor"
  )
  ofit <- suppressWarnings(npreg(
    bws = obw,
    txdat = ox,
    tydat = dat$y,
    exdat = oex,
    gradients = TRUE,
    se = TRUE
  ))
  ooracle <- npreg(
    bws = obw,
    txdat = ox,
    tydat = dat$y,
    exdat = oframes$upper,
    gradients = FALSE
  )$mean - npreg(
    bws = obw,
    txdat = ox,
    tydat = dat$y,
    exdat = oframes$lower,
    gradients = FALSE
  )$mean

  expect_equal(ofit$grad[, 2L], ooracle, tolerance = 1e-10)
  expect_true(all(is.finite(ofit$gerr[, 2L])))
  expect_true(all(ofit$gerr[, 2L] >= 0))
  expect_true(any(ofit$gerr[, 2L] > 0))
})

test_that("general-LP conditional and quantile effects use target endpoints", {
  dat <- a1_effect_data()
  yframe <- data.frame(y = dat$y)
  ey <- data.frame(y = dat$y[seq_len(nrow(dat$ex))])
  common <- list(
    xdat = dat$x,
    ydat = yframe,
    bws = c(0.38, 0.42, 0.20),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 2L,
    basis = "tensor"
  )

  dbw <- do.call(npcdensbw, common)
  dens <- suppressWarnings(npcdens(
    bws = dbw,
    txdat = dat$x,
    tydat = yframe,
    exdat = dat$ex,
    eydat = ey,
    gradients = TRUE
  ))
  dens.oracle <- npcdens(
    bws = dbw,
    txdat = dat$x,
    tydat = yframe,
    exdat = dat$frames$upper,
    eydat = ey,
    gradients = FALSE
  )$condens - npcdens(
    bws = dbw,
    txdat = dat$x,
    tydat = yframe,
    exdat = dat$frames$lower,
    eydat = ey,
    gradients = FALSE
  )$condens
  expect_equal(dens$congrad[, 2L], dens.oracle, tolerance = 1e-10)
  expect_true(all(is.na(dens$congerr[, 2L])))

  fbw <- do.call(npcdistbw, common)
  dist <- suppressWarnings(npcdist(
    bws = fbw,
    txdat = dat$x,
    tydat = yframe,
    exdat = dat$ex,
    eydat = ey,
    gradients = TRUE
  ))
  dist.oracle <- npcdist(
    bws = fbw,
    txdat = dat$x,
    tydat = yframe,
    exdat = dat$frames$upper,
    eydat = ey,
    gradients = FALSE
  )$condist - npcdist(
    bws = fbw,
    txdat = dat$x,
    tydat = yframe,
    exdat = dat$frames$lower,
    eydat = ey,
    gradients = FALSE
  )$condist
  expect_equal(dist$congrad[, 2L], dist.oracle, tolerance = 1e-10)
  expect_true(all(is.na(dist$congerr[, 2L])))

  qfit <- suppressWarnings(npqreg(
    bws = fbw,
    txdat = dat$x,
    tydat = yframe,
    exdat = dat$ex,
    tau = 0.45,
    gradients = TRUE
  ))
  qoracle <- npqreg(
    bws = fbw,
    txdat = dat$x,
    tydat = yframe,
    exdat = dat$frames$upper,
    tau = 0.45,
    gradients = FALSE
  )$quantile - npqreg(
    bws = fbw,
    txdat = dat$x,
    tydat = yframe,
    exdat = dat$frames$lower,
    tau = 0.45,
    gradients = FALSE
  )$quantile
  expect_equal(
    as.numeric(qfit$quantgrad[, 2L]),
    as.numeric(qoracle),
    tolerance = 1e-10
  )
  expect_true(all(is.na(qfit$quantgerr[, 2L])))
})

test_that("general-LP conmode and nplsqreg effects use public target endpoints", {
  dat <- a1_effect_data()
  cls <- factor(
    ifelse(
      dat$y + 0.25 * dat$x$z > stats::median(dat$y),
      "high",
      "low"
    ),
    levels = c("low", "high")
  )
  clsframe <- data.frame(cls = cls)
  cmbw <- npcdensbw(
    xdat = dat$x,
    ydat = clsframe,
    bws = c(0.20, 0.42, 0.20),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 2L,
    basis = "tensor"
  )
  mode.fit <- suppressWarnings(npconmode(
    bws = cmbw,
    txdat = dat$x,
    tydat = clsframe,
    exdat = dat$ex,
    probabilities = TRUE,
    gradients = TRUE,
    level = "high",
    proper = TRUE
  ))
  mode.upper <- npconmode(
    bws = cmbw,
    txdat = dat$x,
    tydat = clsframe,
    exdat = dat$frames$upper,
    probabilities = TRUE,
    gradients = FALSE,
    level = "high",
    proper = TRUE
  )
  mode.lower <- npconmode(
    bws = cmbw,
    txdat = dat$x,
    tydat = clsframe,
    exdat = dat$frames$lower,
    probabilities = TRUE,
    gradients = FALSE,
    level = "high",
    proper = TRUE
  )
  mode.oracle <- mode.upper$probabilities[, "high"] -
    mode.lower$probabilities[, "high"]

  expect_equal(
    mode.fit$probability.gradients[, 2L],
    mode.oracle,
    tolerance = 1e-10
  )
  expect_equal(
    rowSums(mode.fit$probabilities),
    rep(1, nrow(dat$ex)),
    tolerance = 1e-12
  )

  ls.fit <- suppressWarnings(nplsqreg(
    txdat = dat$x,
    tydat = dat$y,
    exdat = dat$ex,
    tau = 0.5,
    bws = c(0.42, 0.20),
    bandwidth.compute = FALSE,
    scale = rep.int(1, nrow(dat$x)),
    regtype = "lp",
    degree = 2L,
    basis = "tensor",
    gradients = TRUE,
    se = TRUE
  ))
  ls.oracle <- nplsqreg(
    bws = ls.fit$bws,
    txdat = dat$x,
    tydat = dat$y,
    exdat = dat$frames$upper,
    gradients = FALSE,
    se = FALSE
  )$quantile - nplsqreg(
    bws = ls.fit$bws,
    txdat = dat$x,
    tydat = dat$y,
    exdat = dat$frames$lower,
    gradients = FALSE,
    se = FALSE
  )$quantile

  expect_equal(ls.fit$quantgrad[, 2L], ls.oracle, tolerance = 1e-10)
  expect_true(all(is.finite(ls.fit$quantgerr[, 2L])))
  expect_true(all(ls.fit$quantgerr[, 2L] >= 0))
  expect_true(any(ls.fit$quantgerr[, 2L] > 0))
})
