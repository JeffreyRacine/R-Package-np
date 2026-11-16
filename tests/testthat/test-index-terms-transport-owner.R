test_that('index transport changes only ordinary private terms ownership', {
  owner <- environment(get('.np_index_dispatch', asNamespace('npRmpi')))
  wire <- get('.np_index_terms_wire', owner)
  role <- get('.np_index_role_wire', owner)
  restore <- get('.np_index_bandwidth_terms_restore', owner)
  e <- new.env(parent = baseenv()); e$padding <- raw(800000L)
  f <- ~x + g; environment(f) <- e; tt <- terms(f)
  p <- wire(tt)
  expect_identical(environment(p), emptyenv())
  expect_identical(environment(tt), e)
  expect_identical(deparse(p), deparse(tt))
  for (light in list(baseenv(), globalenv(), emptyenv())) {
    t <- tt; environment(t) <- light; expect_identical(wire(t), t)
  }
  custom <- tt; class(custom) <- c('custom_terms', class(tt))
  expect_identical(wire(custom), custom)
  x <- data.frame(x = 1:4, g = factor(c('a','b','a','b')))
  attr(x, '.np.index.terms') <- tt
  z <- role(x)
  expect_identical(z$x, x$x); expect_identical(z$g, x$g)
  expect_identical(environment(attr(x, '.np.index.terms')), e)
  b <- list(index.design = list(terms = p), .np.native.training = list(xdat = z),
            call = as.call(list(as.name('npindexbw'), xdat = z)))
  r <- restore(b, tt)
  expect_identical(r$index.design$terms, tt)
  expect_identical(attr(r$.np.native.training$xdat, '.np.index.terms'), tt)
  expect_identical(attr(r$call$xdat, '.np.index.terms'), tt)
  other <- terms(~different); b$index.design$terms <- other
  expect_identical(restore(b, tt)$index.design$terms, other)
  expect_identical(restore(b, NULL), b)
  expect_lt(length(serialize(z, NULL)), length(serialize(x, NULL)) - 790000L)
})

test_that('leased index fits retain the original design terms', {
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves(1L), 'MPI runtime not available')
  on.exit(close_mpi_slaves(), add = TRUE)
  e <- new.env(parent = baseenv()); i <- seq_len(24L)
  e$d <- data.frame(y=cos(i*.3), x=sin(i*.7), g=factor(rep(c('a','b'),12)))
  f <- y ~ x + g; environment(f) <- e
  b <- npindexbw(f, data=e$d, bws=c(1,.5,.4), bandwidth.compute=FALSE)
  expect_identical(environment(b$index.design$terms), e)
  fit <- npindex(b, se=FALSE)
  expect_identical(environment(fit$bws$index.design$terms), e)
  new <- npindex(b, newdata=e$d[1:4,], se=FALSE)
  expect_identical(environment(new$bws$index.design$terms), e)
  expect_equal(predict(fit,newdata=e$d[1:4,]),fitted(new),tolerance=0)
  expect_identical(environment(b$formula), e)
})
