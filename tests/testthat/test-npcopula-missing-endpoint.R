test_that("copula duplicate missing endpoints need no finite rounding scale", {
  f<-getFromNamespace(".npcopula_display_values","npRmpi")
  g<-list(row.order=1:4,display.rows=c(1L,3L),display.cell=c(1L,1L,2L,2L))
  eps<-.Machine$double.eps
  for(missing in c(NA_real_,NaN,Inf,-Inf)) {
    v<-c(missing,missing,1,1+eps);scale<-c(abs(missing),abs(missing),2,2)
    expect_identical(f(v,g,roundoff.scale=scale),v[c(1L,3L)])
    exact<-c(missing,missing,1,1)
    expect_identical(f(exact,g,exact=TRUE,roundoff.scale=scale),exact[c(1L,3L)])
    perm<-c(4L,2L,3L,1L)
    gp<-list(row.order=1:4,display.rows=c(4L,3L),display.cell=g$display.cell[perm])
    expect_identical(f(v[perm],gp,roundoff.scale=scale[perm]),v[c(1L,3L)])
    v[4L]<-1+128*eps
    expect_error(f(v,g,roundoff.scale=scale),"conflicting",fixed=TRUE)
  }
  for(v in list(c(NA_real_,NaN,1,1),c(Inf,-Inf,1,1),c(NA_real_,1,1,1)))
    expect_error(f(v,g,roundoff.scale=abs(v)),"conflicting",fixed=TRUE)
  for(scale in list(c(1,NA_real_,1,1),c(1,Inf,1,1),c(1,-1,1,1),1:2,"x"))
    expect_error(f(rep(1,4),g,roundoff.scale=scale),"operand scale",fixed=TRUE)
  expect_error(f(c(NA,NA,1,1),g,roundoff.scale=c(NA,NA,1,NA)),"operand scale",fixed=TRUE)
})
