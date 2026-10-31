test_that("hat factor payloads use retained response levels without recoding numbers", {
  ns <- asNamespace("np")
  untangle <- get("untangle", ns)
  response <- get(".np_hat_response", ns)
  for (labels in list(c("10", "30", "50", "70"), c("a", "b", "c", "d"))) {
    original <- ordered(labels[1:3], levels = labels)
    dati <- untangle(data.frame(y = original))
    payload <- factor(labels[1:3], levels = rev(labels[1:3]))
    expected <- if (labels[1] == "10") c(10, 30, 50) else c(1, 2, 3)
    expect_identical(as.numeric(response(payload, dati)), expected)
    expect_error(response(factor("unknown"), dati), "unknown factors", fixed = TRUE)
  }
  value <- matrix(1:6, 3, 2)
  expect_identical(response(value, NULL), value)
  payload <- factor(c("a", "b"))
  expect_identical(response(payload, untangle(data.frame(y = 1:2))), payload)
})

test_that("explicit and cached hat factor application matches its numeric operator", {

  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = seq(.05, .95, length.out = 24))
  y <- ordered(rep(c(10, 30, 50), 8))
  b <- npplregbw(xdat = data.frame(cat = y), ydat = seq_along(y),
                  zdat = x, bws = matrix(.3, 2, 1),
                  bandwidth.compute = FALSE)$bw[[2L]]
  H <- npreghat(b, txdat = x)
  expected <- as.vector(H %*% rep(c(10, 30, 50), 8))
  expect_equal(as.vector(npreghat(b, txdat = x, y = y, output = "apply")),
               expected, tolerance = 2e-11)
  expect_equal(predict(H, y = y, output = "apply"), expected, tolerance = 2e-11)
})


test_that("single-index primary responses require an explicit numeric scale", {
  withr::local_options(np.messages = FALSE)
  x <- data.frame(x = seq(.05, .95, length.out = 24), z = sin(seq_len(24)))
  y <- ordered(rep(c("a", "b", "c"), 8))
  expect_error(npindexbw(xdat=x, ydat=y, bws=c(1,.4,.3),
                         bandwidth.compute=FALSE), "must be a numeric")
  scores <- rep(c(10, 30, 50), 8)
  b <- npindexbw(xdat=x, ydat=scores, bws=c(1,.4,.3), bandwidth.compute=FALSE)
  expect_identical(fitted(npindex(b,se=FALSE)),
                   fitted(npindex(b,txdat=x,tydat=scores,se=FALSE)))
})
