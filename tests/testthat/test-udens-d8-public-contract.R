test_that("D8 default selection preserves explicit controls and other owners", {
  b <- structure(list(type="fixed", ckertype="gaussian", ckerorder=2L,
    ckerbound="none", ndim=2L, ncon=2L), class="bandwidth")
  pick <- function(control=NULL, center="bias-corrected", method="asymptotic", alpha=.05)
    .np_udens_d8_default(b, control, center, method, alpha)
  expect_true(pick())
  expect_false(pick(np_udens_bias_control()))
  expect_false(pick(method="bootstrap"))
  expect_false(pick(method="none"))
  expect_false(pick(center="estimate"))
  for (a in c(1e-5,1e-4,.5,NA_real_)) expect_false(pick(alpha=a))
  for (q in 1:4) { b$ndim <- b$ncon <- q; expect_identical(pick(),q<=3L) }
  b$ndim <- 2L; b$ncon <- 1L; expect_false(pick())
})

test_that("D8 row layout is exact and both workspaces retire together", {
  s <- .np_udens_d8_store(5L)
  on.exit(.np_udens_d8_store_close(s),add=TRUE)
  for(j in 1:3) .np_udens_d8_store_append(s, as.double(1:5+j))
  reference <- .np_udens_d8_store_rows(s,2L,4L)
  .np_udens_d8_store_for_draws(s,2L)
  expect_identical(.np_udens_d8_store_rows(s,2L,4L),reference)
  expect_error(.np_udens_d8_store_for_draws(s,2L),"transition")
  paths <- c(s$path,s$row.path)
  .np_udens_d8_store_close(s)
  expect_false(any(file.exists(paths)))
})

test_that("D8 native boundary rejects invalid matrix and column inputs", {
  x <- matrix(c(-1,0,1),3L);g<-matrix(c(-.5,.5),2L)
  call <- function(x, j=1L) .Call(C_np_plot_udens_d8_column,x,g,.5,c(.2,.2),j)
  expect_error(call(matrix(1:3,3L)),"double matrices")
  expect_error(call(x,0L),"dimensions or column")
  expect_error(call(x,NA_integer_),"dimensions or column")
  expect_error(call(x,1.0),"dimensions or column")
  x[1L] <- Inf
  expect_error(call(x),"nonfinite")
})


test_that("D8 completed-rank layout preserves admission and explicit tiles", {
  state <- list(n=100L,q=2L,max.bytes=256*1024^2)
  for(processes in c(1L,2L,4L)) {
    budget <- .np_udens_d8_pilot_budget(state,2500L,9999L,3L,NULL,32L,256L,processes)
    actual <- .np_udens_d8_pilot_draw_budget(state,2500L,417L,budget)
    expect_lte(actual$bytes,state$max.bytes)
    expect_gte(actual$tile,32L)
    expect_lte(actual$tile,2500L)
    for(delta in c(-1,0,1)) {
      cut <- state;cut$max.bytes <- budget$bytes+delta
      if(delta<0) expect_error(.np_udens_d8_pilot_budget(cut,2500L,9999L,3L,NULL,32L,256L,processes),"byte budget")
      else expect_no_error(.np_udens_d8_pilot_budget(cut,2500L,9999L,3L,NULL,32L,256L,processes))
    }
    explicit <- .np_udens_d8_pilot_budget(state,2500L,9999L,3L,32L,32L,256L,processes)
    expect_identical(.np_udens_d8_pilot_draw_budget(state,2500L,417L,explicit)$tile,32L)
    expect_error(.np_udens_d8_pilot_draw_budget(state,2500L,0L,budget),"factor rank")
  }
})
