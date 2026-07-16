.glp_complete_package <- "np"

.glp_complete_nsfun <- function(name) {
  get(name, envir = asNamespace(.glp_complete_package), inherits = FALSE)
}

.glp_complete_local <- function(expr) {
  eval(substitute(expr), envir = parent.frame())
}

test_that("complete GLP terms and dimensions are pinned", {
  build <- .glp_complete_nsfun("npBuildLpTerms")
  ncol_basis <- .glp_complete_nsfun("npLpBasisNcol")
  dim_bs <- .glp_complete_nsfun("dimBS")

  fixtures <- list(
    `2,1` = list(
      degree = c(2L, 1L),
      terms = rbind(
        c(0L, 0L), c(1L, 0L), c(2L, 0L), c(0L, 1L), c(1L, 1L)
      )
    ),
    `2,2` = list(
      degree = c(2L, 2L),
      terms = rbind(
        c(0L, 0L), c(1L, 0L), c(2L, 0L),
        c(0L, 1L), c(1L, 1L), c(0L, 2L)
      )
    ),
    `2,0` = list(
      degree = c(2L, 0L),
      terms = rbind(c(0L, 0L), c(1L, 0L), c(2L, 0L))
    ),
    `1,2` = list(
      degree = c(1L, 2L),
      terms = rbind(
        c(0L, 0L), c(1L, 0L), c(0L, 1L), c(1L, 1L), c(0L, 2L)
      )
    )
  )

  for (fixture in fixtures) {
    degree <- fixture$degree
    expected <- fixture$terms
    expect_identical(unname(build(degree, basis = "glp")), expected)
    expect_equal(as.numeric(ncol_basis("glp", degree)), nrow(expected),
                 tolerance = 0)
    expect_equal(
      as.numeric(dim_bs("glp", kernel = TRUE, degree = degree,
                        segments = rep.int(1L, length(degree)))),
      nrow(expected) - 1L,
      tolerance = 0
    )
  }

  expect_equal(nrow(build(c(2L, 1L), basis = "additive")), 4L,
               tolerance = 0)
  expect_equal(nrow(build(c(2L, 1L), basis = "tensor")), 6L,
               tolerance = 0)
})

test_that("complete GLP representation label names shifted Legendre", {
  label <- .glp_complete_nsfun("npLpBasisRepresentationLabel")
  expect_identical(label(TRUE, "glp"),
                   "Shifted Legendre (degree-graded orthonormal)")
  expect_identical(label(TRUE, "additive"), "Bernstein")
  expect_identical(label(FALSE, "glp"), "Raw")
})

test_that("complete GLP shifted-Legendre levels and derivatives are analytic", {
  mypoly <- .glp_complete_nsfun("mypoly")
  w_lp <- .glp_complete_nsfun("W.lp")

  train <- c(-2, -1, 0, 1, 3)
  eval <- c(-1.5, 0.5, 2)
  xrange <- diff(range(train))
  z <- 2 * (eval - min(train)) / xrange - 1
  dz <- 2 / xrange

  expected0 <- cbind(
    sqrt(3) * z,
    sqrt(5) * (3 * z^2 - 1) / 2,
    sqrt(7) * (5 * z^3 - 3 * z) / 2
  )
  expected1 <- cbind(
    rep(sqrt(3) * dz, length(z)),
    sqrt(5) * 3 * z * dz,
    sqrt(7) * (15 * z^2 - 3) * dz / 2
  )
  expected2 <- cbind(
    rep(0, length(z)),
    rep(sqrt(5) * 3 * dz^2, length(z)),
    sqrt(7) * 15 * z * dz^2
  )
  expected3 <- cbind(
    rep(0, length(z)),
    rep(0, length(z)),
    rep(sqrt(7) * 15 * dz^3, length(z))
  )

  actual <- lapply(0:3, function(order) {
    mypoly(
      x = train,
      ex = eval,
      degree = 3L,
      gradient.compute = order > 0L,
      r = order,
      Bernstein = TRUE,
      complete.glp = TRUE
    )
  })

  expect_equal(actual[[1L]], expected0, tolerance = 5e-14)
  expect_equal(actual[[2L]], expected1, tolerance = 5e-14)
  expect_equal(actual[[3L]], expected2, tolerance = 5e-14)
  expect_equal(actual[[4L]], expected3, tolerance = 5e-14)
  expect_equal(
    mypoly(x = train, ex = eval, degree = 3L, gradient.compute = TRUE,
           r = 4L, Bernstein = TRUE, complete.glp = TRUE),
    matrix(0, nrow = length(eval), ncol = 3L),
    tolerance = 0
  )

  train2 <- data.frame(x1 = train, x2 = seq(1, 5, length.out = length(train)))
  eval2 <- data.frame(x1 = eval, x2 = c(1.5, 3, 4.5))
  d1 <- w_lp(train2, eval2, degree = c(2L, 1L),
             gradient.vec = c(1L, 0L), basis = "glp",
             bernstein.basis = TRUE)
  d2 <- w_lp(train2, eval2, degree = c(2L, 1L),
             gradient.vec = c(0L, 1L), basis = "glp",
             bernstein.basis = TRUE)

  z1 <- 2 * (eval2$x1 - min(train2$x1)) / diff(range(train2$x1)) - 1
  z2 <- 2 * (eval2$x2 - min(train2$x2)) / diff(range(train2$x2)) - 1
  phi1.x1 <- sqrt(3) * z1
  phi1.x2 <- sqrt(3) * z2
  dphi1.x1 <- sqrt(3) * 2 / diff(range(train2$x1))
  dphi1.x2 <- sqrt(3) * 2 / diff(range(train2$x2))

  expect_equal(unname(d1[, "1.1"]), dphi1.x1 * phi1.x2,
               tolerance = 5e-14)
  expect_equal(unname(d2[, "1.1"]), phi1.x1 * dphi1.x2,
               tolerance = 5e-14)
})

test_that("raw and shifted-Legendre GLP fits agree on public functionals", {
  set.seed(20260716)
  n <- 64L
  tx <- data.frame(
    x1 = runif(n, -0.8, 0.8),
    x2 = runif(n, -0.7, 0.9)
  )
  y <- 0.4 + 0.7 * tx$x1 - 0.3 * tx$x1^2 + 0.5 * tx$x2 +
    0.2 * tx$x2^2 + 0.8 * tx$x1 * tx$x2 +
    0.03 * sin(seq_len(n) * 1.7)
  ex <- data.frame(
    x1 = seq(-0.6, 0.6, length.out = 9L),
    x2 = seq(-0.5, 0.7, length.out = 9L)
  )

  make_bw <- function(bernstein, degree) {
    .glp_complete_local(npregbw(
      xdat = tx,
      ydat = y,
      regtype = "lp",
      degree = degree,
      degree.select = "manual",
      basis = "glp",
      bernstein.basis = bernstein,
      bws = c(0.45, 0.45),
      bandwidth.compute = FALSE
    ))
  }
  fit_with <- function(bw) {
    .glp_complete_local(npreg(
      txdat = tx,
      tydat = y,
      exdat = ex,
      bws = bw,
      gradients = TRUE,
      gradient.order = c(2L, 1L)
    ))
  }

  for (degree in list(c(2L, 1L), c(2L, 2L))) {
    raw <- fit_with(make_bw(FALSE, degree))
    shifted <- fit_with(make_bw(TRUE, degree))

    expect_equal(raw$mean, shifted$mean, tolerance = 1e-11,
                 info = paste(degree, collapse = ","))
    expect_equal(raw$grad, shifted$grad, tolerance = 1e-10,
                 info = paste(degree, collapse = ","))
    expect_equal(raw$merr, shifted$merr, tolerance = 1e-11,
                 info = paste(degree, collapse = ","))
    expect_equal(raw$gerr, shifted$gerr, tolerance = 1e-10,
                 info = paste(degree, collapse = ","))
  }
})
