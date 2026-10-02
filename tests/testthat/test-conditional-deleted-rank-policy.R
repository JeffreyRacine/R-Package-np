# Structural protection complements separately retained native boundary probes.
test_that("conditional deleted-rank policy is common and degree neutral", {
  roots <- unique(c(test_path("..", ".."), Sys.getenv("R_PACKAGE_SOURCE"),
                    Sys.getenv("R_PACKAGE_DIR"), getwd()))
  roots <- roots[nzchar(roots)]
  roots <- roots[file.exists(file.path(roots, "src", "conditional_rank_admission.h"))]
  skip_if(!length(roots), "package C sources unavailable")
  read <- function(name) paste(readLines(file.path(roots[1], "src", name),
                                         warn = FALSE), collapse = "\n")
  policy <- read("conditional_rank_admission.h")
  global <- read("conditional_global_qr.h")
  source <- read("jksum.c")
  expect_match(policy, "np_conditional_cold_rank(workspace, p)", fixed = TRUE)
  expect_match(policy, "rank_upper_bound = 0;", fixed = TRUE)
  expect_match(policy, "np_lp_solve_workspace_solve_adjoint_ranked", fixed = TRUE)
  expect_false(grepl("degree", policy, fixed = TRUE))
  # The adopted local-design certificate precedes the original-coordinate
  # deficient-rank solve; the obsolete raw-Gram certificate is not required.
  expect_match(global, "np_cqr_local_row(", fixed = TRUE)
  expect_match(global, "if(rank==NP_CQR_LOCAL_FULL)", fixed = TRUE)
  expect_match(global, "if(rank!=NP_CQR_LOCAL_DEFICIENT)", fixed = TRUE)
  expect_match(global, "np_conditional_solve_adjoint_ranked(&policy,k,1,1./n,0,&d)",
               fixed = TRUE)
  expect_match(global, "np_conditional_solve_adjoint_ranked(", fixed = TRUE)
  expect_match(source, "np_conditional_solve_adjoint_ranked(", fixed = TRUE)
})
