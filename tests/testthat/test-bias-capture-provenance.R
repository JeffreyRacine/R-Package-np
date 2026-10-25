test_that("regression capture compares model state independently of call provenance", {
  same <- np:::.np_plot_regression_capture_same_bws
  a <- list(bw=c(.5,.3),degree=c(1L,2L),type="fixed",call=quote(npregbw(y~x)))
  attr(a$call,".Environment") <- environment()
  b <- a;attr(b$call,".Environment") <- NULL
  expect_true(same(a,b))
  before <- a;expect_true(same(a,a));expect_identical(a,before)
  for(field in c("bw","degree","type")) {
    changed <- b;changed[[field]] <- switch(field,bw=c(.5,.4),degree=c(1L,3L),type="adaptive_nn")
    expect_false(same(a,changed))
  }
  b$new.computational.field <- TRUE
  expect_false(same(a,b))
  expect_false(same(NULL,a))
})
