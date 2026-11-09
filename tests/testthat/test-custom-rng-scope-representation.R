test_that("scope representation checks refuse uncertifiable custom generators", {
  set.seed(123);entry<-.Random.seed
  builtin<-c(10403L,rep(0L,625L))
  expect_null(.np_rng_scope_check(builtin))
  uniform<-builtin;uniform[1L]<-10405L
  expect_error(.np_rng_scope_check(uniform),"user-supplied generator state is not fully represented",fixed=TRUE)
  expect_error(.np_rng_scope_check(uniform,uniform.only=TRUE),"user-supplied generator state is not fully represented",fixed=TRUE)
  normal<-builtin;normal[1L]<-10303L
  expect_error(.np_rng_scope_check(normal),"user-supplied generator state is not fully represented",fixed=TRUE)
  expect_null(.np_rng_scope_check(normal,uniform.only=TRUE))
  expect_identical(.Random.seed,entry)
})
