test_that("conditional objects retain properization controls for consumers", {
  old <- options(np.messages=FALSE)
  on.exit(options(old),add=TRUE)
  set.seed(920318)
  d <- data.frame(x=runif(30,-1,1))
  d$y <- sin(4*d$x)+rnorm(30,sd=.2)
  e <- data.frame(x=c(-.4,.1,.6),y=c(-.3,.2,.4))
  for (family in c("npcdens","npcdist")) for (formula in c(FALSE,TRUE)) {
    args <- if (formula) list(formula=y~x,data=d) else
      list(xdat=d["x"],ydat=d["y"])
    bw <- do.call(get(paste0(family,"bw")),c(args,list(bws=c(.3,.3),
      bandwidth.compute=FALSE,regtype="ll")))
    ctrl <- list(mode="slice",slice.grid.size=41L,apply="both")
    fit <- do.call(get(family),list(bws=bw,newdata=e,proper=TRUE,
                                     proper.control=ctrl))
    expect_true(fit$proper.applied)
    expect_identical(fit$proper.control$mode,"slice")
    expect_identical(fit$proper.control$slice.grid.size,41L)
    expect_equal(predict(fit,newdata=e),fitted(fit),tolerance=1e-12)
    expect_equal(predict(fit,newdata=e,proper.control=ctrl),fitted(fit),
                 tolerance=1e-12)
    raw <- do.call(get(family),list(bws=bw,newdata=e,proper=FALSE))
    expect_equal(predict(fit,newdata=e,proper=FALSE),fitted(raw),tolerance=0)
    expect_error(predict(fit,newdata=e,proper.control=list()),
                 "repeated|grid|slice")
    legacy <- fit
    legacy$proper.control <- NULL
    expect_error(predict(legacy,newdata=e),"repeated|grid|slice")
    # Retention contains only normalized scalar policy, never data or a grid.
    expect_true(all(lengths(fit$proper.control)==1L))
    expect_lt(as.numeric(object.size(fit$proper.control)),4096)
    replay <- .np_conditional_replay_proper(fit,list())
    expect_identical(replay$proper.control,fit$proper.control)
    strict <- .np_conditional_replay_proper(fit,list(),prediction=TRUE)
    expect_true(strict$proper.control$fail.on.unsupported)
    explicit <- .np_conditional_replay_proper(fit,
      list(proper.control=list(fail.on.unsupported=FALSE)),prediction=TRUE)
    expect_false(explicit$proper.control$fail.on.unsupported)
  }
})
