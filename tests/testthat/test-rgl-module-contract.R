.rgl_test_package <- if (isNamespaceLoaded("npRmpi")) "npRmpi" else "np"

test_that("canonical rgl controls preserve established parsing and precedence", {
  controls <- getFromNamespace(".np_plot_rgl_controls", .rgl_test_package)
  collect <- getFromNamespace(".np_plot_collect_rgl_args", .rgl_test_package)
  extract <- getFromNamespace(".np_plot_extract_prefixed_args", .rgl_test_package)
  merge.legend <- getFromNamespace(".np_plot_merge_rgl_legend_control", .rgl_test_package)

  dots <- list(
    fov = 70,
    zoom = 0.8,
    rgl.view3d.fov = 45,
    rgl.par3d.windowRect = c(1, 2, 301, 302),
    rgl.grid3d.col = "grey70",
    rgl.widget.width = 720,
    rgl.legend3d.x = "bottomleft",
    rgl.surface3d.alpha = 0.3,
    rgl.points3d.size = 4
  )

  result <- controls(dots, legend = list(show = TRUE, inset = 0.1))

  expect_identical(
    result$persp3d,
    collect(dots, "rgl.persp3d", "rgl.persp3d.")
  )
  expect_identical(
    result$view3d,
    collect(dots, "rgl.view3d", "rgl.view3d.")
  )
  expect_identical(result$view3d$fov, 45)
  expect_identical(result$view3d$zoom, 0.8)
  expect_identical(
    result$par3d,
    collect(dots, "rgl.par3d", "rgl.par3d.")
  )
  expect_identical(
    result$grid3d,
    collect(dots, "rgl.grid3d", "rgl.grid3d.")
  )
  expect_identical(
    result$widget,
    collect(dots, "rgl.widget", "rgl.widget.")
  )
  expect_identical(
    result$legend3d,
    merge.legend(
      collect(dots, "rgl.legend3d", "rgl.legend3d."),
      list(show = TRUE, inset = 0.1)
    )
  )
  expect_identical(
    result$surface3d,
    extract(dots, "rgl.surface3d.")
  )
  expect_identical(
    result$points3d,
    extract(dots, "rgl.points3d.")
  )
})

test_that("canonical rgl controls retain every legend condition contract", {
  controls <- getFromNamespace(".np_plot_rgl_controls", .rgl_test_package)

  expect_identical(controls(list(), TRUE)$legend3d, list())
  expect_identical(controls(list(), FALSE)$legend3d, list(plot = FALSE))
  expect_identical(controls(list(), NA)$legend3d, list(plot = FALSE))
  expect_identical(controls(list(), NULL)$legend3d, list(plot = FALSE))
  expect_identical(controls(list(), "topleft")$legend3d, list(x = "topleft"))
  expect_identical(
    controls(list(), list(show = TRUE, x = "bottomright"))$legend3d,
    list(x = "bottomright")
  )
  expect_error(controls(list(), logical(0L)), "legend must be")
  expect_error(controls(list(), list(show = 1)), "legend\\$show")
})

test_that("supported plot routes use one rgl control owner", {
  root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  module <- file.path(root, "R", "np.plot.rgl.R")
  skip_if_not(file.exists(module), "source R files unavailable in installed test context")

  module.source <- paste(readLines(module, warn = FALSE), collapse = "\n")
  expect_match(
    module.source,
    ".np_plot_rgl_controls <- function",
    fixed = TRUE
  )

  engine.files <- Sys.glob(file.path(root, "R", "np.plot.engine*.R"))
  engine.source <- paste(unlist(lapply(engine.files, readLines, warn = FALSE)),
                         collapse = "\n")
  expect_false(grepl("rgl.persp3d.user.args", engine.source, fixed = TRUE))
  expect_false(grepl(".np_plot_collect_rgl_args", engine.source, fixed = TRUE))
  expect_gte(
    lengths(regmatches(engine.source,
                       gregexpr(".np_plot_rgl_controls", engine.source,
                                fixed = TRUE))),
    7L
  )

  direct.source <- paste(
    readLines(file.path(root, "R", "np.copula.R"), warn = FALSE),
    readLines(file.path(root, "R", "np.plot.methods.R"), warn = FALSE),
    collapse = "\n"
  )
  expect_match(direct.source, ".np_plot_rgl_controls", fixed = TRUE)
})
