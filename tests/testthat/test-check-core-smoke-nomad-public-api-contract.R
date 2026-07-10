nomad_public_methods <- list(
  npcdensbw = c("conbandwidth", "default"),
  npcdistbw = c("condbandwidth", "default"),
  npindexbw = "default",
  nplsqregbw = "default",
  npregbw = c("default", "rbandwidth"),
  npscoefbw = "default",
  npudensbw = c("bandwidth", "default"),
  npudistbw = c("dbandwidth", "default")
)

nomad_method_list <- function() {
  out <- list()
  for (generic in names(nomad_public_methods)) {
    for (class in nomad_public_methods[[generic]]) {
      out[[paste(generic, class, sep = ".")]] <- getS3method(generic, class)
    }
  }
  out
}

stub_from <- function(fun) {
  as.function(c(formals(fun), quote(match.call(expand.dots = FALSE))),
              envir = baseenv())
}

test_that("NOMAD-capable bandwidth methods expose keyword-only nomad.opts", {
  for (method in nomad_method_list()) {
    fml <- names(formals(method))
    expect_true("..." %in% fml)
    expect_true("nomad.opts" %in% fml)
    expect_gt(match("nomad.opts", fml), match("...", fml))
  }
})

test_that("nomad.opts argument matching is exact and non-positional", {
  for (method in nomad_method_list()) {
    stub <- stub_from(method)

    exact <- do.call(stub, list(nomad.opts = list(MAX_BB_EVAL = "5")))
    expect_true("nomad.opts" %in% names(exact))
    expect_false("nomad.opts" %in% names(exact[["..."]]))

    partial <- do.call(stub, list(nomad.opt = list(MAX_BB_EVAL = "5")))
    expect_false("nomad.opts" %in% names(partial))
    expect_true("nomad.opt" %in% names(partial[["..."]]))

    positional <- do.call(stub, as.list(seq_len(length(formals(stub)) + 2L)))
    expect_false("nomad.opts" %in% names(positional))
    expect_true("..." %in% names(positional))
  }
})

test_that("nomad.opts normalization rejects non-list controls before solver entry", {
  normalize <- getFromNamespace(".np_nomad_normalize_user_opts", "np")

  expect_identical(normalize(NULL, "npregbw"), list())
  expect_identical(normalize(list(MAX_BB_EVAL = "5"), "npregbw"),
                   list(MAX_BB_EVAL = "5"))
  expect_error(normalize(1, "npregbw"), "'nomad.opts' must be a list", fixed = TRUE)

  for (method in nomad_method_list()) {
    expect_error(
      method(nomad.opts = 1),
      "'nomad.opts' must be a list",
      fixed = TRUE
    )
  }
})
