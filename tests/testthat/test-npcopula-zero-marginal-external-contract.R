.st6a_capture_warnings <- function(expr) {
  messages <- character()
  value <- withCallingHandlers(
    force(expr),
    warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = messages)
}

.st6a_copula_fixture <- function() {
  z <- c(seq(0, 0.01, length.out = 10L),
         seq(0.99, 1, length.out = 10L))
  dat <- data.frame(
    x = z,
    y = z + seq(-1e-4, 1e-4, length.out = length(z))
  )
  h <- c(0.01, 0.01)
  list(
    dat = dat,
    h = h,
    dbw = npudensbw(
      dat = dat, ckertype = "epanechnikov", ckerorder = 2L,
      bwtype = "fixed", bandwidth.compute = FALSE, bws = h
    ),
    fbw = npudistbw(
      dat = dat, ckertype = "epanechnikov", ckerorder = 2L,
      bwtype = "fixed", bandwidth.compute = FALSE, bws = h
    )
  )
}

.st6a_copula_facts <- function(object, fixture) {
  out <- as.data.frame(object)
  xgrid <- out[, c("x", "y"), drop = FALSE]
  marginal <- lapply(seq_len(2L), function(j) {
    name <- c("x", "y")[[j]]
    bw <- npudensbw(
      dat = fixture$dat[name],
      ckertype = "epanechnikov", ckerorder = 2L,
      bwtype = "fixed", bandwidth.compute = FALSE,
      bws = fixture$h[j]
    )
    fitted(npudens(
      bws = bw, tdat = fixture$dat[name], edat = xgrid[name]
    ))
  })
  joint <- fitted(npudens(
    bws = fixture$dbw, tdat = fixture$dat, edat = xgrid
  ))
  data.frame(
    x = xgrid$x, y = xgrid$y, joint = joint,
    mx = marginal[[1L]], my = marginal[[2L]],
    zero.marginal = marginal[[1L]] == 0 | marginal[[2L]] == 0,
    copula = out$copula
  )
}

test_that("external copula-density zero marginals are rowwise undefined", {
  if (exists("spawn_mpi_slaves", mode = "function")) {
    if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }

  fixture <- .st6a_copula_fixture()
  supported <- .st6a_capture_warnings(npcopula(
    bws = fixture$dbw, data = fixture$dat,
    u = data.frame(x = 0.49, y = 0.49),
    n.quasi.inv = 200L, er.quasi.inv = 0
  ))
  expect_length(supported$warnings, 0L)
  expect_equal(fitted(supported$value), 9.2085385703235,
               tolerance = 1e-14)

  grid <- .st6a_capture_warnings(npcopula(
    bws = fixture$dbw, data = fixture$dat,
    u = data.frame(
      x = c(0.49, 0.5, 0.5001),
      y = c(0.49, 0.5, 0.5001)
    ),
    n.quasi.inv = 200L, er.quasi.inv = 0
  ))
  facts <- .st6a_copula_facts(grid$value, fixture)
  expect_length(grid$warnings, 1L)
  expect_match(grid$warnings, "^\\[np(Rmpi)?\\]")
  expect_match(grid$warnings, "5 of 9 external evaluation rows",
               fixed = TRUE)
  expect_match(grid$warnings, "copula density ratio is undefined",
               fixed = TRUE)
  expect_identical(is.na(facts$copula), facts$zero.marginal)
  expect_true(any(facts$joint == 0 & !facts$zero.marginal))
  expect_true(all(facts$copula[facts$joint == 0 &
                               !facts$zero.marginal] == 0))

  all.undefined <- .st6a_capture_warnings(npcopula(
    bws = fixture$dbw, data = fixture$dat,
    u = data.frame(x = rep(0.5, 4L), y = rep(0.5, 4L)),
    n.quasi.inv = 200L, er.quasi.inv = 0
  ))
  expect_true(all(is.na(fitted(all.undefined$value))))
  expect_length(all.undefined$warnings, 1L)
  expect_match(all.undefined$warnings, "16 of 16 external evaluation rows",
               fixed = TRUE)
  expect_match(all.undefined$warnings, "and 8 others", fixed = TRUE)

  predicted <- .st6a_capture_warnings(predict(
    supported$value,
    u = data.frame(x = 0.5, y = 0.5),
    output = "object"
  ))
  expect_length(predicted$warnings, 1L)
  expect_s3_class(predicted$value, "npcopula")
  expect_true(is.na(fitted(predicted$value)))
  expect_true(is.na(as.data.frame(predicted$value)$copula))

  distribution <- .st6a_capture_warnings(npcopula(
    bws = fixture$fbw, data = fixture$dat,
    u = data.frame(x = 0.5, y = 0.5),
    n.quasi.inv = 200L, er.quasi.inv = 0
  ))
  expect_length(distribution$warnings, 0L)
  expect_identical(as.data.frame(distribution$value)$copula, 0.5)
})
