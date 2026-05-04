library(npRmpi)
options(np.messages = FALSE)

.np_demo_src <- Sys.getenv("NP_DEMO_SRC", "")
.np_demo_family <- c(if (nzchar(.np_demo_src)) file.path(.np_demo_src, "..", "inst", "demo_family_nptests.R"),
                     system.file("demo_family_nptests.R", package = "npRmpi"))
.np_demo_family <- .np_demo_family[nzchar(.np_demo_family) & file.exists(.np_demo_family)]
source(.np_demo_family[[1L]])
nptest_demo_source_utils()

nslaves <- np_demo_n(default = 1L, floor = 1L,
                     exact_env = "NP_DEMO_NSLAVES",
                     frac_env = "NP_DEMO_NSLAVES_FRAC")
npRmpi.init(nslaves = nslaves)

nptest_demo_run_matrix("npdeneqtest", "session")
