test_that("npregbw exposes nomad.opts explicitly without moving existing npRmpi formals", {
  methods <- list(
    default = getS3method("npregbw", "default"),
    rbandwidth = getS3method("npregbw", "rbandwidth")
  )

  for (method in methods) {
    fml <- names(formals(method))
    expect_true("..." %in% fml)
    expect_true("nomad.opts" %in% fml)
    expect_lt(match("nomad.opts", fml), match("...", fml))
  }
})

test_that("npregbw preserves exact named nomad.opts binding in npRmpi", {
  stub_from <- function(fun) {
    as.function(c(formals(fun), quote(match.call(expand.dots = FALSE))),
                envir = baseenv())
  }

  stub <- stub_from(getS3method("npregbw", "default"))

  exact <- do.call(stub, list(xdat = 1, ydat = 2, nomad.opts = list(MAX_BB_EVAL = "5")))
  expect_true("nomad.opts" %in% names(exact))
  expect_false("nomad.opts" %in% names(exact[["..."]]))

  positional <- do.call(stub, as.list(seq_len(match("nomad.opts", names(formals(stub))))))
  expect_true("nomad.opts" %in% names(positional))
})

test_that("npregbw normalizes only list-like NOMAD option controls in npRmpi", {
  normalize <- getFromNamespace(".np_nomad_normalize_user_opts", "npRmpi")

  expect_identical(normalize(NULL, "npregbw"), list())
  expect_identical(normalize(list(MAX_BB_EVAL = "5"), "npregbw"),
                   list(MAX_BB_EVAL = "5"))
  expect_error(normalize(1, "npregbw"), "'nomad.opts' must be a list", fixed = TRUE)

  default_method <- getS3method("npregbw", "default")
  expect_error(
    default_method(
      xdat = data.frame(x = seq_len(5)),
      ydat = seq_len(5),
      bandwidth.compute = FALSE,
      nomad.opts = 1
    ),
    "'nomad.opts' must be a list",
    fixed = TRUE
  )
})
