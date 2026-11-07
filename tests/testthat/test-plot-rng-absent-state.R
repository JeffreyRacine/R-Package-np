test_that("absent custom uniform seed is not mistaken for complete state", {
  predicate <- get(".np_plot_replay_uniform_state_absent",asNamespace("npRmpi"))
  present <- exists(".Random.seed",.GlobalEnv,inherits=FALSE)
  prior <- if(present)get(".Random.seed",.GlobalEnv,inherits=FALSE) else NULL
  clear <- function() if(exists(".Random.seed",.GlobalEnv,inherits=FALSE))
    rm(list=".Random.seed",envir=.GlobalEnv)
  on.exit({clear();if(present)assign(".Random.seed",prior,.GlobalEnv)},add=TRUE)
  custom <- c("user-supplied","Inversion","Rejection")
  builtin <- c("Mersenne-Twister","Inversion","Rejection")
  clear()
  expect_true(predicate(custom))
  expect_false(predicate(builtin))
  expect_false(exists(".Random.seed",.GlobalEnv,inherits=FALSE))
  assign(".Random.seed",10405L,.GlobalEnv)
  expect_true(predicate(custom))
  expect_false(predicate(builtin))
  full <- c(10405L,1L,2L);assign(".Random.seed",full,.GlobalEnv)
  expect_false(predicate(custom));expect_identical(.Random.seed,full)
  reads <- 0L;clear()
  makeActiveBinding(".Random.seed",function(value){reads <<- reads+1L;stop("seed binding evaluated")},.GlobalEnv)
  expect_false(predicate(custom));expect_false(predicate(builtin))
  expect_identical(reads,0L)
})

test_that("both preflight owners give absent and header seeds the same decision", {
  ns <- asNamespace("npRmpi")
  ordinary <- get(".np_plot_rng_replay_preflight",ns)
  factor <- get(".np_plot_replay_factor_preflight",ns)
  e <- new.env(parent=environment(ordinary))
  e$RNGkind <- function()c("user-supplied","Inversion","Rejection")
  environment(ordinary) <- e;environment(factor) <- e
  e$.np_plot_replay_factor_preflight <- factor
  present <- exists(".Random.seed",.GlobalEnv,inherits=FALSE)
  prior <- if(present)get(".Random.seed",.GlobalEnv,inherits=FALSE)else NULL
  clear <- function()if(exists(".Random.seed",.GlobalEnv,inherits=FALSE))rm(list=".Random.seed",envir=.GlobalEnv)
  on.exit({clear();if(present)assign(".Random.seed",prior,.GlobalEnv)},add=TRUE)
  for(mode in c("absent","header","full")) {
    clear()
    if(mode!="absent")assign(".Random.seed",if(mode=="header")10405L else c(10405L,1L,2L),.GlobalEnv)
    for(family in c("regression","cdf","density")) {
      b <- structure(list(type="fixed",ncon=1L,nuno=0L,nord=0L,
        ckertype="gaussian",ckerorder=2L,ckerbound="none"),
        class=switch(family,regression="rbandwidth",cdf="dbandwidth",density="bandwidth"))
      method <- get(switch(family,regression=".np_plot_rbandwidth_engine",
        cdf=".np_plot_dbandwidth_engine",density=".np_plot_bandwidth_engine"),ns)
      d <- list(plot.errors.center="bias-corrected",plot.errors.method="asymptotic",
        plot.errors.type="simultaneous",neval=if(family=="density")513L else 33L)
      if(family!="regression")d$bias.control <- if(family=="cdf")
        np_udist_bias_control("cdf") else np_udens_bias_control("joint-scale")
      if(mode=="full" && family=="regression")expect_silent(ordinary(method,b,d)) else
        if(mode=="full")expect_error(ordinary(method,b,d),"user-supplied generator state cannot be certified",fixed=TRUE) else
        expect_error(ordinary(method,b,d),"custom uniform state absent",fixed=TRUE)
      if(family!="regression") {
        if(mode=="full")expect_error(factor(method,b,d),"user-supplied generator state cannot be certified",fixed=TRUE) else
          expect_error(factor(method,b,d),"custom uniform state absent",fixed=TRUE)
        d$neval <- if(family=="density")512L else 32L
        expect_silent(factor(method,b,d))
      }
    }
    if(mode=="absent")expect_false(exists(".Random.seed",.GlobalEnv,inherits=FALSE))
  }
})
