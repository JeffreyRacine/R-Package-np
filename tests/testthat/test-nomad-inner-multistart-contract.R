test_that("nomad inner multistart detects named dots and accepts active NOMAD search", {
  out <- .np_nomad_validate_inner_multistart(
    call_names = c("xdat", "ydat", "..."),
    dot.args = list(nomad.nmulti = 2L),
    nomad.nmulti = 0L,
    regtype = "lp",
    automatic.degree.search = TRUE,
    search.engine = "nomad+powell"
  )

  expect_true(isTRUE(out$named))
  expect_identical(out$nmulti, 2L)
})

test_that("nomad inner multistart fails fast when named dots are outside active NOMAD search", {
  expect_error(
    .np_nomad_validate_inner_multistart(
      call_names = c("xdat", "ydat", "..."),
      dot.args = list(nomad.nmulti = 1L),
      nomad.nmulti = 0L,
      regtype = "lc",
      automatic.degree.search = FALSE,
      search.engine = "nomad+powell"
    ),
    "nomad.nmulti is only supported"
  )
})
