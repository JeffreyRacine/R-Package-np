test_that(".docall resolves function names from caller frame", {
  plus1 <- function(x) x + 1L
  expect_identical(.docall("plus1", list(2L)), 3L)
})

test_that(".docall falls back to .GlobalEnv function lookup", {
  nm <- ".__npRmpi_docall_global_fun"
  had <- exists(nm, envir = .GlobalEnv, inherits = FALSE)
  if (had) {
    old <- get(nm, envir = .GlobalEnv, inherits = FALSE)
  }
  assign(nm, function(x) x + 2L, envir = .GlobalEnv)
  on.exit({
    if (had) {
      assign(nm, old, envir = .GlobalEnv)
    } else if (exists(nm, envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = nm, envir = .GlobalEnv)
    }
  }, add = TRUE)

  expect_identical(.docall(nm, list(1L)), 3L)
})

test_that(".docall lookup path uses caller-frame get0 first", {
  fn.body <- paste(deparse(body(.docall), width.cutoff = 500L), collapse = " ")
  expect_match(fn.body, "get0\\(fname, envir = envir, mode = \"function\", inherits = TRUE\\)")
  expect_match(fn.body, "get0\\(fname, envir = \\.GlobalEnv, mode = \"function\", inherits = FALSE\\)")
})
