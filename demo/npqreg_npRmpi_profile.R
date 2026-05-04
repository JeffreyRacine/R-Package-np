mpi.bcast.cmd(np.mpi.initialize(), caller.execute = TRUE)
mpi.bcast.cmd(options(np.messages = FALSE), caller.execute = TRUE)

.np_demo_src <- Sys.getenv("NP_DEMO_SRC", "")
.np_demo_family <- c(if (nzchar(.np_demo_src)) file.path(.np_demo_src, "..", "inst", "demo_family_npaux.R"),
                     system.file("demo_family_npaux.R", package = "npRmpi"))
.np_demo_family <- .np_demo_family[nzchar(.np_demo_family) & file.exists(.np_demo_family)]
source(.np_demo_family[[1L]])
mpi.bcast.Robj2slave(.np_demo_family)
mpi.bcast.cmd(source(.np_demo_family[[1L]]), caller.execute = FALSE)
npaux_demo_source_utils()
npaux_demo_run_matrix("npqreg", "profile")

mpi.bcast.cmd(mpi.quit(), caller.execute = TRUE)
