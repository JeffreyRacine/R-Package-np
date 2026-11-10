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


test_that("scope header validation respects binomial-kind capability", {
  entry <- if (exists('.Random.seed', .GlobalEnv, inherits=FALSE)) .Random.seed else NULL
  for (count in c(3L, 4L)) {
    check <- .np_rng_scope_check
    environment(check) <- list2env(list(.np_rng_kind_count=function() count),
      parent=environment(check))
    old <- c(10403L, rep.int(1L, 625L)); modern <- old; modern[1L] <- 110403L
    expect_null(check(old))
    if (count == 4L) expect_null(check(modern)) else
      expect_error(check(modern), 'unrecognized seed representation')
    for (header in c(-1L, NA_integer_, 210403L, 120403L, 110803L)) {
      bad <- old; bad[1L] <- header
      expect_error(check(bad), 'seed|representation')
    }
    expect_error(check(old[-2L]), 'unrecognized seed representation')
    for (prefix in c(0L, 100000L)) {
      bad <- old; bad[1L] <- prefix + 10203L
      expect_error(check(bad), 'Box-Muller')
      bad[1L] <- prefix + 10405L
      expect_error(check(bad, uniform.only=TRUE), 'user-supplied')
      bad[1L] <- prefix + 10303L
      expect_error(check(bad), 'user-supplied')
    }
    old[20L] <- NA_integer_
    expect_null(check(old))
  }
  expect_identical(if (exists('.Random.seed', .GlobalEnv, inherits=FALSE)) .Random.seed else NULL, entry)
})

test_that("RNG capability inspection leaves an absent seed absent", {
  had <- exists('.Random.seed', .GlobalEnv, inherits=FALSE)
  old <- if(had) .Random.seed else NULL
  on.exit(if(had) assign('.Random.seed',old,.GlobalEnv),add=TRUE)
  if(had) rm('.Random.seed',envir=.GlobalEnv)
  expect_identical(.np_rng_kind_count(), length(formals(base::RNGkind)))
  expect_false(exists('.Random.seed',.GlobalEnv,inherits=FALSE))
})
