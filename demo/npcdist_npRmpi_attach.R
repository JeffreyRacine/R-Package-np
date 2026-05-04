library(npRmpi)
npRmpi.init(mode = "attach", comm = 1, autodispatch = TRUE)
options(np.messages = FALSE)

if (mpi.comm.rank(0L) == 0L) {
  .np_demo_src <- Sys.getenv("NP_DEMO_SRC", "")
  .np_demo_family <- c(if (nzchar(.np_demo_src)) file.path(.np_demo_src, "..", "inst", "demo_family_npcdist.R"),
                       system.file("demo_family_npcdist.R", package = "npRmpi"))
  .np_demo_family <- .np_demo_family[nzchar(.np_demo_family) & file.exists(.np_demo_family)]
  source(.np_demo_family[[1L]])

  npcdist_demo_run_matrix("attach")
  npRmpi.quit(mode = "attach", comm = 1)
}
