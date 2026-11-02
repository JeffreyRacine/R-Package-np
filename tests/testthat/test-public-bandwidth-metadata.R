test_that('metadata wire keeps R syntax, attributes and canonical identity', {
  wire <- getFromNamespace('.npRmpi_bw_metadata_wire','npRmpi')
  fp <- getFromNamespace('.npRmpi_autodispatch_fingerprint','npRmpi')
  make <- function(bytes) {
    unrelated <- raw(bytes)
    f <- y ~ x + z
    d <- data.frame(y=1:5,x=6:10,z=11:15)
    frame <- model.frame(f,d)
    call <- substitute(npregbw(formula=FORMULA,...=DOTS),
      list(FORMULA=f,DOTS=pairlist(bws=.2,bandwidth.compute=FALSE)))
    environment(call) <- environment()
    list(call=call,formula=f,frame=frame,terms=attr(frame,'terms'),
      extra=pairlist(a=NULL,b=quote(expr=)))
  }
  original <- make(0L); projected <- wire(original)
  expect_identical(fp(original),fp(projected))
  expect_type(projected$call[['...']],'pairlist')
  expect_type(projected$extra,'pairlist')
  expect_identical(names(projected$extra),c('a','b'))
  expect_identical(names(attributes(original$terms)),names(attributes(projected$terms)))
  expect_identical(environment(projected$formula),emptyenv())
  expect_false(identical(environment(original$formula),emptyenv()))
  expect_identical(length(serialize(wire(make(1000000L)),NULL)),
                   length(serialize(projected,NULL)))
  expect_identical(wire(projected),projected)
})

test_that('transport substitution distinguishes formulas from ordinary calls', {
  replace <- getFromNamespace('.npRmpi_autodispatch_replace_tmps','npRmpi')
  owner <- new.env(parent=baseenv());owner$a<-1L
  f <- y ~ .__npRmpi_autod_x_1
  environment(f)<-owner
  expected <- quote(f(c(1,2),NULL))
  plain <- quote(f(.__npRmpi_autod_x_1,NULL))
  got <- replace(plain,list(.__npRmpi_autod_x_1=quote(c(1,2))))
  expect_identical(got,expected)
  expect_identical(replace(f,list(.__npRmpi_autod_x_1=quote(c(1,2)))),f)
  wrapped <- as.call(list(as.name('f'),f,NULL))
  expect_identical(replace(wrapped,list(.__npRmpi_autod_x_1=quote(c(1,2)))),wrapped)
})

test_that('vector tau metadata preserves scalar calls and finishes the aggregate', {
  apply <- getFromNamespace('.npRmpi_bw_formula_apply','npRmpi')
  f <- y ~ x
  child <- structure(list(call=quote(nplsqregbw(xdat=x,ydat=y,tau=.3)),
    rows.omit=NA_integer_,nobs.omit=0L,formula=NULL,ynames='y'),class='lsqregressionbandwidth')
  aggregate <- child;aggregate$terms<-terms(f);aggregate$tau.bws<-list(child)
  plan <- list(producer='nplsqregbw',formula=f,terms=terms(f),ynames='y',
    omit=structure(2L,class='exclude'),call=quote(nplsqregbw(y~x,tau=c(.3,.7))),
    scalar.children=TRUE)
  got <- apply(aggregate,plan)
  expect_identical(got$call,plan$call)
  expect_identical(tail(names(got),1L),'terms')
  expect_identical(got$tau.bws[[1L]]$terms,terms(f))
  expect_identical(got$tau.bws[[1L]]$call,child$call)
  expect_identical(got$rows.omit,2L)
  expect_identical(got$tau.bws[[1L]]$rows.omit,NA_integer_)
  plan$scalar.child<-TRUE
  got <- apply(child,plan)
  expect_identical(got$call,child$call)
  expect_identical(got$rows.omit,child$rows.omit)
  expect_identical(got$formula,f)
})
