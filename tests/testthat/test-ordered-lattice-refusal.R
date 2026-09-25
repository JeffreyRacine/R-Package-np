test_that("raw LR cumulative refuses ambiguity without rounding fractions", {
  check <- .np_ordered_lr_cumulative_contract
  expect_error(check(c(1023.1,1025.1,1026.1)), "numerically ambiguous")
  for(s in list(c(.5,1,2), c(23.1,25.1,26.1), c(1023.5,1025.5),
                1e14+c(0,1.125,2.125), c(1023.1,1025.1,1026.6)))
    expect_true(check(s))
  expect_error(.np_ordered_levels_contract(c(1023.1,1024.1), "nliracine"),
               "cannot be established reliably")
  expect_true(.np_ordered_levels_contract(c(.5,1,2), "liracine"))
  expect_true(.np_ordered_levels_contract(c(.5,1,2), "racineliyan"))
})
