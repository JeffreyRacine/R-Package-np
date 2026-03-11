library(testthat)
library(npRmpi)

local({
  on.exit({
    try(npRmpi.quit(force = TRUE), silent = TRUE)
    try(mpi.finalize(), silent = TRUE)
    if ("package:npRmpi" %in% search()) {
      try(detach("package:npRmpi", unload = TRUE, character.only = TRUE), silent = TRUE)
    }
    if ("npRmpi" %in% loadedNamespaces()) {
      try(unloadNamespace("npRmpi"), silent = TRUE)
    }
  }, add = TRUE)
  test_check("npRmpi")
})
