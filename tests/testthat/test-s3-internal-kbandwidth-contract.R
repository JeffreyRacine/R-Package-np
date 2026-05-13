test_that("internal kbandwidth S3 methods are registered but not exported", {
  ns <- asNamespace(getNamespaceName(environment(npksum)))
  s3 <- getNamespaceInfo(ns, "S3methods")

  expect_false("kbandwidth" %in% getNamespaceExports(ns))
  expect_true(any(s3[, 1L] == "kbandwidth" & s3[, 2L] == "numeric"))
  expect_true(any(s3[, 1L] == "kbandwidth" & s3[, 2L] == "integer"))
  expect_true(any(s3[, 1L] == "kbandwidth" & s3[, 2L] == "default"))
  expect_true(any(s3[, 1L] == "print" & s3[, 2L] == "kbandwidth"))
  expect_true(any(s3[, 1L] == "as.double" & s3[, 2L] == "kbandwidth"))

  expect_true(!is.null(getS3method("print", "kbandwidth", optional = TRUE)))
  expect_true(!is.null(getS3method("as.double", "kbandwidth", optional = TRUE)))
})

test_that("internal kbandwidth objects print and coerce consistently", {
  kbandwidth <- getFromNamespace("kbandwidth", getNamespaceName(environment(npksum)))
  untangle <- getFromNamespace("untangle", getNamespaceName(environment(npksum)))

  x <- data.frame(x = seq(0, 1, length.out = 5))
  kbw <- kbandwidth(
    bw = 0.5,
    xdati = untangle(x),
    xnames = names(x),
    nobs = nrow(x)
  )

  expect_s3_class(kbw, "kbandwidth")
  expect_equal(as.double(kbw), 0.5)
  out <- capture.output(print(kbw))
  expect_true(any(grepl("Kernel Sum Bandwidth", out, fixed = TRUE)))
  expect_true(any(grepl("Bandwidth Type: Fixed", out, fixed = TRUE)))
})
