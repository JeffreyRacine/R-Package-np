test_that("public bootstrap-count owners expose B with the historical defaults", {
  owners <- list(
    npcmstest = npcmstest,
    npqcmstest = npqcmstest,
    npdeneqtest = npdeneqtest,
    npdeptest = npdeptest,
    npsdeptest = npsdeptest,
    npsigtest = npsigtest,
    npsigtest.rbandwidth = getS3method("npsigtest", "rbandwidth"),
    npindex = npindex,
    npindex.sibandwidth = getS3method("npindex", "sibandwidth"),
    npsymtest = npsymtest,
    npunitest = npunitest
  )

  for (owner in owners) {
    owner.formals <- formals(owner)
    expect_true("B" %in% names(owner.formals))
    expect_false("boot.num" %in% names(owner.formals))
    expect_identical(owner.formals$B, 399)
  }

  expect_identical(names(formals(npsigtest)), c("bws", "...", "B"))
  expect_identical(names(formals(npindex)), c("bws", "...", "B"))
  expected.positions <- c(
    npcmstest = 9L, npqcmstest = 11L, npdeneqtest = 5L,
    npdeptest = 5L, npsdeptest = 5L, npsigtest = 3L,
    npsigtest.rbandwidth = 4L, npindex = 3L,
    npindex.sibandwidth = 6L, npsymtest = 3L, npunitest = 5L
  )
  observed.positions <- vapply(
    owners,
    function(owner) match("B", names(formals(owner))),
    integer(1L)
  )
  expect_identical(observed.positions, expected.positions)
})

test_that("retired boot.num is not accepted by public bootstrap-count owners", {
  expect_error(npcmstest(boot.num = 9), "use 'B' instead", fixed = TRUE)
  expect_error(npqcmstest(boot.num = 9), "use 'B' instead", fixed = TRUE)
  expect_error(npdeneqtest(boot.num = 9), "use 'B' instead", fixed = TRUE)
  expect_error(npdeptest(boot.num = 9), "unused argument", fixed = TRUE)
  expect_error(npsdeptest(boot.num = 9), "unused argument", fixed = TRUE)
  expect_error(npsigtest(boot.num = 9), "use 'B' instead", fixed = TRUE)
  expect_error(npindex(boot.num = 9), "use 'B' instead", fixed = TRUE)
  expect_error(
    getS3method("npindex", "sibandwidth")(boot.num = 9),
    "use 'B' instead", fixed = TRUE
  )
  expect_error(npsymtest(boot.num = 9), "use 'B' instead", fixed = TRUE)
  expect_error(npunitest(boot.num = 9), "use 'B' instead", fixed = TRUE)
})
