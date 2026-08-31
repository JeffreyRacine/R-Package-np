test_that("partial-linear bootstrap helpers materialize zero-ridge-first policy", {
  helper.names <- c(
    ".np_plreg_weighted_coef",
    ".np_inid_boot_from_plreg_exact",
    ".np_inid_boot_from_plreg_frozen",
    ".np_inid_boot_from_plreg"
  )
  defaults <- vapply(helper.names, function(name) {
    eval(formals(getFromNamespace(name, "np"))$ridge)
  }, numeric(1L))
  expect_identical(unname(defaults), rep(0, length(helper.names)))

  ridge.sequence <- getFromNamespace("npRidgeSequenceFromBase", "np")
  expect_identical(
    ridge.sequence(n.train = 6L, ridge.base = 0, cap = 1),
    seq.int(from = 0, to = 1, by = 1 / 6)
  )
})
