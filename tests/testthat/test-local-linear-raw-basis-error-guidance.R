library(npRmpi)

local_linear_guidance_dens_msg <- "C_np_density_conditional_bw: optimizer returned a fixed-bandwidth candidate with invalid raw objective"
local_linear_guidance_dist_msg <- "C_np_distribution_conditional_bw: optimizer returned a fixed-bandwidth candidate with invalid raw objective"

local_linear_guidance_ll_spec <- function(ncon = 1L) {
  npRmpi:::npCanonicalConditionalRegSpec(regtype = "ll", ncon = ncon)
}

local_linear_guidance_lp_spec <- function(ncon = 1L) {
  npRmpi:::npCanonicalConditionalRegSpec(
    regtype = "lp",
    degree = rep.int(1L, ncon),
    bernstein.basis = FALSE,
    ncon = ncon
  )
}

test_that("raw degree-1 request helper matches legacy ll and explicit raw lp", {
  expect_true(npRmpi:::npIsRawDegreeOneConditionalRequest("ll"))
  expect_true(npRmpi:::npIsRawDegreeOneConditionalRequest(
    "lp",
    degree = 1L,
    bernstein.basis = FALSE
  ))
  expect_true(npRmpi:::npIsRawDegreeOneConditionalRequest(
    "lp",
    degree = c(1L, 1L),
    bernstein.basis = FALSE
  ))
  expect_false(npRmpi:::npIsRawDegreeOneConditionalRequest(
    "lp",
    degree = 1L,
    bernstein.basis = TRUE
  ))
  expect_false(npRmpi:::npIsRawDegreeOneConditionalRequest(
    "lp",
    degree = 2L,
    bernstein.basis = FALSE
  ))
  expect_false(npRmpi:::npIsRawDegreeOneConditionalRequest("lc"))
})

test_that("local-linear raw-basis guidance is targeted to conditional density ll cv.ls", {
  expect_error(
    npRmpi:::npWithLocalLinearRawBasisSearchError(
      stop(local_linear_guidance_dens_msg, call. = FALSE),
      where = "npcdensbw",
      spec = local_linear_guidance_ll_spec(),
      bwmethod = "cv.ls",
      ncon = 1L
    ),
    regexp = "npcdensbw\\(\\) local-linear cv\\.ls failed while using the canonical raw degree-1 basis"
  )

  expect_error(
    npRmpi:::npWithLocalLinearRawBasisSearchError(
      stop(local_linear_guidance_dens_msg, call. = FALSE),
      where = "npcdensbw",
      spec = local_linear_guidance_lp_spec(),
      bwmethod = "cv.ls",
      ncon = 1L
    ),
    regexp = "npcdensbw\\(\\) local-linear cv\\.ls failed while using the canonical raw degree-1 basis"
  )
})

test_that("local-linear raw-basis guidance is targeted to conditional distribution ll cv.ls", {
  expect_error(
    npRmpi:::npWithLocalLinearRawBasisSearchError(
      stop(local_linear_guidance_dist_msg, call. = FALSE),
      where = "npcdistbw",
      spec = local_linear_guidance_ll_spec(),
      bwmethod = "cv.ls",
      ncon = 1L
    ),
    regexp = "npcdistbw\\(\\) local-linear cv\\.ls failed while using the canonical raw degree-1 basis"
  )

  expect_error(
    npRmpi:::npWithLocalLinearRawBasisSearchError(
      stop(local_linear_guidance_dist_msg, call. = FALSE),
      where = "npcdistbw",
      spec = local_linear_guidance_lp_spec(),
      bwmethod = "cv.ls",
      ncon = 1L
    ),
    regexp = "npcdistbw\\(\\) local-linear cv\\.ls failed while using the canonical raw degree-1 basis"
  )
})

test_that("local-linear raw-basis guidance leaves non-targeted errors unchanged", {
  expect_error(
    npRmpi:::npWithLocalLinearRawBasisSearchError(
      stop(local_linear_guidance_dist_msg, call. = FALSE),
      where = "npcdensbw",
      spec = local_linear_guidance_ll_spec(),
      bwmethod = "cv.ls",
      ncon = 1L
    ),
    regexp = local_linear_guidance_dist_msg
  )

  expect_error(
    npRmpi:::npWithLocalLinearRawBasisSearchError(
      stop(local_linear_guidance_dens_msg, call. = FALSE),
      where = "npcdensbw",
      spec = local_linear_guidance_ll_spec(),
      bwmethod = "cv.ml",
      ncon = 1L
    ),
    regexp = local_linear_guidance_dens_msg
  )

  expect_error(
    npRmpi:::npWithLocalLinearRawBasisSearchError(
      stop(local_linear_guidance_dens_msg, call. = FALSE),
      where = "npcdensbw",
      spec = local_linear_guidance_ll_spec(ncon = 0L),
      bwmethod = "cv.ls",
      ncon = 0L
    ),
    regexp = local_linear_guidance_dens_msg
  )

  bernstein_spec <- npRmpi:::npCanonicalConditionalRegSpec(
    regtype = "lp",
    degree = 1L,
    bernstein.basis = TRUE,
    ncon = 1L
  )
  expect_error(
    npRmpi:::npWithLocalLinearRawBasisSearchError(
      stop(local_linear_guidance_dens_msg, call. = FALSE),
      where = "npcdensbw",
      spec = bernstein_spec,
      bwmethod = "cv.ls",
      ncon = 1L
    ),
    regexp = local_linear_guidance_dens_msg
  )
})
