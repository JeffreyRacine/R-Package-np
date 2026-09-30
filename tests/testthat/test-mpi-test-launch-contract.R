test_that("test launch opt-in preserves MPI ownership and independent bootstrap", {
  launcher <- Sys.which("mpiexec")
  skip_if(.Platform$OS.type!="unix" || !nzchar(launcher),
          "explicit Unix MPI launcher unavailable")
  spec <- npRmpi_test_launch_spec
  expect_identical(spec("Rscript","body",launcher="",inherited="PMI_FD"),
                   list(command="Rscript",args="body"))
  child <- spec("/path/Rscript","body",launcher,
                inherited=c("PMI_FD","PMI_HOSTNAME"))
  expect_identical(child$command,"/usr/bin/env")
  expect_identical(child$args,
    c("-u","PMI_FD","-u","PMI_HOSTNAME",shQuote(launcher),"-n","1",
      shQuote("/path/Rscript"),"body"))
  for(command in c(launcher,"/path/R","/bin/echo"))
    expect_identical(spec(command,c("-n","3"),launcher,inherited=character())$args,
                     c(shQuote(command),"-n","3"))
  for(key in c("PMI_UNKNOWN","PMIX_RANK","OMPI_COMM_WORLD_RANK"))
    expect_error(spec("Rscript","body",launcher,inherited=key),
                 "unqualified inherited MPI bootstrap",fixed=TRUE)
  expect_error(spec("Rscript","body",launcher,inherited=character(),
                    env="PMI_FD=6"),"cannot override MPI bootstrap",fixed=TRUE)
})
