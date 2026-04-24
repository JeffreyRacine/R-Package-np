library(np)

gen_bw_sel_str <- getFromNamespace("genBwSelStr", "np")

test_that("bandwidth summaries map base-run multistart labels to one", {
  expect_match(
    gen_bw_sel_str(list(fval = 1.23, ifval = 0, num.feval = NA, num.feval.fast = NA)),
    "achieved on multistart 1",
    fixed = TRUE
  )

  expect_match(
    gen_bw_sel_str(list(fval = 1.23, ifval = 2, num.feval = NA, num.feval.fast = NA)),
    "achieved on multistart 2",
    fixed = TRUE
  )
})

test_that("bandwidth summaries ignore malformed ifval and prefer NOMAD best restart", {
  expect_false(grepl(
    "achieved on multistart",
    gen_bw_sel_str(list(fval = 1.23, ifval = 4.81432545643515e-06, num.feval = NA, num.feval.fast = NA)),
    fixed = TRUE
  ))

  expect_match(
    gen_bw_sel_str(list(
      fval = 1.23,
      ifval = 4.81432545643515e-06,
      nomad.best.restart = 3L,
      num.feval = NA,
      num.feval.fast = NA
    )),
    "achieved on multistart 3",
    fixed = TRUE
  )
})

test_that("single-start progress avoids legacy redraw output", {
  set.seed(1)
  x <- runif(20)
  y <- x + rnorm(20, sd = 0.1)

  out <- capture.output(
    npcdensbw(y ~ x, regtype = "ll", bwtype = "adaptive_nn", nmulti = 1, itmax = 1)
  )

  expect_false(any(grepl("Multistart 1 of 1", out, fixed = TRUE)))
  expect_false(any(grepl("Multistart 1 of 0", out, fixed = TRUE)))
})
