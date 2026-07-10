test_that("npregbw exposes nomad.opts as an explicit keyword-only public control", {
  methods <- list(
    default = getS3method("npregbw", "default"),
    rbandwidth = getS3method("npregbw", "rbandwidth")
  )

  for (method in methods) {
    fml <- names(formals(method))
    expect_true("..." %in% fml)
    expect_true("nomad.opts" %in% fml)
    expect_gt(match("nomad.opts", fml), match("...", fml))
  }
})

test_that("npregbw nomad.opts argument matching is exact and non-positional", {
  stub_from <- function(fun) {
    as.function(c(formals(fun), quote(match.call(expand.dots = FALSE))),
                envir = baseenv())
  }

  stub <- stub_from(getS3method("npregbw", "default"))

  exact <- do.call(stub, list(xdat = 1, ydat = 2, nomad.opts = list(MAX_BB_EVAL = "5")))
  expect_true("nomad.opts" %in% names(exact))
  expect_false("nomad.opts" %in% names(exact[["..."]]))

  partial <- do.call(stub, list(xdat = 1, ydat = 2, nomad.opt = list(MAX_BB_EVAL = "5")))
  expect_false("nomad.opts" %in% names(partial))
  expect_true("nomad.opt" %in% names(partial[["..."]]))

  positional <- do.call(stub, as.list(seq_len(length(formals(stub)) + 2L)))
  expect_false("nomad.opts" %in% names(positional))
  expect_true("..." %in% names(positional))
})

test_that("npregbw normalizes only list-like NOMAD option controls", {
  normalize <- getFromNamespace(".np_nomad_normalize_user_opts", "np")

  expect_identical(normalize(NULL, "npregbw"), list())
  expect_identical(normalize(list(MAX_BB_EVAL = "5"), "npregbw"),
                   list(MAX_BB_EVAL = "5"))
  expect_error(normalize(1, "npregbw"), "'nomad.opts' must be a list", fixed = TRUE)

  expect_error(
    npregbw(
      xdat = data.frame(x = seq_len(5)),
      ydat = seq_len(5),
      bandwidth.compute = FALSE,
      nomad.opts = 1
    ),
    "'nomad.opts' must be a list",
    fixed = TRUE
  )
})
