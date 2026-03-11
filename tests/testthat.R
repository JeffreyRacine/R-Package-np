library(testthat)
library(npRmpi)

local({
  on.exit({
    try(npRmpi.quit(force = TRUE), silent = TRUE)
    try(mpi.finalize(), silent = TRUE)
  }, add = TRUE)
  test_check("npRmpi")
})
