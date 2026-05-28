test_that("fast-path option validators reject silent repair", {
  withr::local_options(np.largeh = TRUE, np.largelambda = FALSE,
                       np.extendednn = TRUE,
                       np.largeh.rel.tol = 0.002,
                       np.disc.upper.rel.tol = 0.02)

  expect_true(np:::npLogicalOption("np.largeh", TRUE))
  expect_false(np:::npLogicalOption("np.largelambda", TRUE))
  expect_true(np:::npExtendedNnEnabled())
  expect_equal(np:::npLargehRelTol(), 0.002)
  expect_equal(np:::npDiscUpperRelTol(), 0.02)

  withr::local_options(np.largeh = 1)
  expect_error(np:::npLogicalOption("np.largeh", TRUE), "TRUE or FALSE")

  withr::local_options(np.largeh = NA)
  expect_error(np:::npLogicalOption("np.largeh", TRUE), "TRUE or FALSE")

  withr::local_options(np.largeh.rel.tol = -1)
  expect_error(np:::npLargehRelTol(), "finite numeric scalar")

  withr::local_options(np.disc.upper.rel.tol = 0.75)
  expect_error(np:::npDiscUpperRelTol(), "finite numeric scalar")
})

test_that("legacy fast-path environment fallbacks do not override R options", {
  old_largeh <- Sys.getenv("NP_LARGEH_REL_TOL", unset = NA_character_)
  old_disc <- Sys.getenv("NP_DISC_UPPER_REL_TOL", unset = NA_character_)
  withr::defer({
    if (is.na(old_largeh)) Sys.unsetenv("NP_LARGEH_REL_TOL")
    else Sys.setenv(NP_LARGEH_REL_TOL = old_largeh)
    if (is.na(old_disc)) Sys.unsetenv("NP_DISC_UPPER_REL_TOL")
    else Sys.setenv(NP_DISC_UPPER_REL_TOL = old_disc)
  })

  Sys.setenv(NP_LARGEH_REL_TOL = "0.05", NP_DISC_UPPER_REL_TOL = "0.25")
  withr::local_options(np.largeh.rel.tol = 0.002,
                       np.disc.upper.rel.tol = 0.02)

  expect_equal(np:::npLargehRelTol(), 0.002)
  expect_equal(np:::npDiscUpperRelTol(), 0.02)
})
