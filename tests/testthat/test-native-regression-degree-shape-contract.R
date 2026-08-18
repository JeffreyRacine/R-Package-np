test_that("native regression objectives reject malformed LP degree vectors", {
  ns <- asNamespace("npRmpi")
  reg.c <- get("npRegtypeToC", ns)(
    regtype = "lp",
    degree = 1L,
    ncon = 1L,
    context = "native degree-shape test"
  )
  myopti <- as.integer(list(
    2L,
    get("IMULTI_FALSE", ns),
    1L,
    get("USE_START_YES", ns),
    get("SF_NORMAL", ns),
    get("BW_FIXED", ns),
    0L,
    get("RE_MIN_FALSE", ns),
    get("IO_MIN_TRUE", ns),
    get("BWM_CVLS", ns),
    get("CKER_GAUSS", ns),
    get("UKER_AIT", ns),
    get("OKER_WANG", ns),
    0L,
    0L,
    1L,
    reg.c$code,
    get("DO_TREE_NO", ns),
    FALSE,
    3L,
    FALSE,
    0L,
    2L,
    TRUE,
    2L
  ))

  expect_error(
    .Call(
      "C_np_regression_bw_eval",
      numeric(),
      numeric(),
      c(0, 1),
      c(0, 1),
      1,
      myopti,
      double(19L),
      0.2,
      1L,
      0L,
      10,
      integer(),
      0L,
      1L,
      -Inf,
      Inf,
      PACKAGE = "npRmpi"
    ),
    "glp_degree length mismatch",
    fixed = TRUE
  )
})
