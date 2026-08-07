test_that("documented rgl prefixes share the canonical renderer allowlist", {
  prefixes <- c("rgl.persp3d", "rgl.view3d", "rgl.par3d", "rgl.grid3d",
                "rgl.widget", "rgl.legend3d")
  expected <- unlist(lapply(prefixes, function(type) {
    paste0(type, ".", npRmpi:::.np_plot_user_arg_names(type))
  }), use.names = FALSE)

  expect_setequal(npRmpi:::.np_plot_rgl_prefixed_arg_names(), expected)
  expect_no_error(npRmpi:::.np_plot_validate_public_dots(
    setNames(as.list(rep(1, length(expected))), expected)
  ))
})

test_that("public rgl prefix validation remains exact and closed", {
  expect_error(
    npRmpi:::.np_plot_validate_public_dots(list(rgl.view3d.not_a_control = 1)),
    "unused plot argument",
    fixed = TRUE
  )
  expect_error(
    npRmpi:::.np_plot_validate_public_dots(list(rgl.view3d. = 1)),
    "unused plot argument",
    fixed = TRUE
  )
  expect_error(
    npRmpi:::.np_plot_validate_public_dots(list(rgl.points3d.size = 4)),
    "unused plot argument",
    fixed = TRUE
  )
})
