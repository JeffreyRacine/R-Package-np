test_that("guide-region streams identify unsupported RNG kinds without changing requests", {
  request <- getFromNamespace(".np_udens_beta_random_request", "npRmpi")
  uniform <- c("Wichmann-Hill", "Marsaglia-Multicarry", "Super-Duper", "Mersenne-Twister", "Knuth-TAOCP", "Knuth-TAOCP-2002", "L'Ecuyer-CMRG")
  normal <- c("Buggy Kinderman-Ramage", "Ahrens-Dieter", "Box-Muller", "Inversion", "Kinderman-Ramage")
  set.seed(71); entry <- .Random.seed
  for (u in uniform) for (n in normal) for (mode in c("normal", "multinomial")) {
    kind <- c(u,n,"Rejection")
    expected <- list(mode=mode,n=10L,B=11L,seed=42L,kind=kind,
                     bytes=110*if(mode=="normal")8 else 4,version=R.version$version.string)
    expect_identical(request(mode,10L,11L,42L,kind),expected)
  }
  expect_error(request("normal",10L,11L,42L,c("user-supplied","Inversion","Rejection")), "beta/Gaussian.*supported built-in RNG")
  expect_error(request("normal",10L,11L,42L,c("Mersenne-Twister","user-supplied","Rejection")), "beta/Gaussian.*supported built-in RNG")
  expect_error(request("normal",10L,0L,42L,c("Mersenne-Twister","Inversion","Rejection")), "invalid beta plot random stream")
  expect_error(request("normal",10L,11L,42L,c("Mersenne-Twister",NA_character_,"Rejection")), "invalid beta plot random stream")
  expect_identical(.Random.seed,entry)
})
