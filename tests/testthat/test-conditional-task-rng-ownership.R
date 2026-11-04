test_that("precomputed conditional task execution preserves statistical advancement", {
  withr::local_seed(613)
  kinds <- RNGkind(); on.exit(do.call(RNGkind, as.list(kinds)), add=TRUE)
  withr::local_options(np.messages=FALSE,np.plot.progress=FALSE)
  for(kind in c("Mersenne-Twister","L'Ecuyer-CMRG")) for(B in c(127L,128L,129L)) {
    RNGkind(kind,"Inversion","Rejection");set.seed(613)
    bank<-matrix(rnorm(3L*B),3L,B);state<-.Random.seed;next.draw<-rnorm(3L)
    set.seed(613);tasks<-.np_conditional_tasks(list(m=3L),B,TRUE)
    expect_identical(do.call(cbind,lapply(tasks,`[[`,"value")),bank)
    expect_identical(attr(tasks,"rng_final_state"),state)
    out<-.np_conditional_dispatch(tasks,list(factor=diag(3)),"gaussian",3L)
    expect_identical(out,bank);expect_identical(.Random.seed,state)
    expect_identical(rnorm(3L),next.draw)
  }
})
