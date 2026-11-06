test_that("bounded Box-Muller regression replay retains actual shared innovations", {
  old.kind <- RNGkind()
  skip_if(identical(old.kind[2L],"Box-Muller"),"test cannot snapshot an incoming hidden normal cache")
  had.seed <- exists(".Random.seed",.GlobalEnv,inherits=FALSE)
  old.seed <- if(had.seed)get(".Random.seed",.GlobalEnv,inherits=FALSE) else NULL
  on.exit({
    do.call(RNGkind,as.list(old.kind))
    if(had.seed)assign(".Random.seed",old.seed,.GlobalEnv) else
      if(exists(".Random.seed",.GlobalEnv,inherits=FALSE))rm(".Random.seed",envir=.GlobalEnv)
  },add=TRUE)
  for(pending in 0:1) {
    reset <- function() {
      RNGkind("Mersenne-Twister","Box-Muller","Rejection");set.seed(71L)
      if(pending)invisible(rnorm(1L))
    }
    reset();z <- matrix(rnorm(7L*33L),7L);end <- .Random.seed;following <- rnorm(5L)
    reset();cache <- npRmpi:::.np_bias_gaussian_cache(7L,33L,"gaussian")
    expect_identical(cache$values,z)
    expect_identical(cache$end,end)
    for(i in 1:2)expect_identical(
      npRmpi:::.np_bias_student_errors(rep(1,7L),33L,"gaussian",cache$values),list(z))
    expect_identical(rnorm(5L),following)
  }
  saved <- .Random.seed
  expect_error(npRmpi:::.np_bias_gaussian_cache(2048L,1025L,"gaussian"),
    "unsafe RNG replay: resolved Gaussian regression",fixed=TRUE)
  expect_error(npRmpi:::.np_bias_gaussian_cache(0L,33L,"gaussian"),
    "invalid Gaussian replay dimensions",fixed=TRUE)
  expect_null(npRmpi:::.np_bias_gaussian_cache(7L,33L,"mammen"))
  expect_identical(.Random.seed,saved)
  RNGkind(normal.kind="Inversion")
  saved <- .Random.seed
  expect_null(npRmpi:::.np_bias_gaussian_cache(2048L,1025L,"gaussian"))
  expect_identical(.Random.seed,saved)
})
