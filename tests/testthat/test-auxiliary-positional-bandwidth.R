test_that("quantile and classification defaults retain positional manual bandwidths", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages=FALSE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(x=seq(-1,1,length.out=18))
  d <- data.frame(x=x$x,y=sin(seq_len(18)),cls=factor(rep(letters[1:3],6)))
  for(family in c("npqreg","npconmode")) {
    ff <- get(family)
    y <- if(family=="npqreg") d$y else d$cls
    bw <- if(family=="npqreg") c(.4,.35) else c(.2,.35)
    named <- ff(bws=bw,txdat=x,tydat=y,se=FALSE)
    for(shape in c("native","all.positional","partial","wrapper")) {
      state <- new.env(parent = emptyenv()); state$n <- 0L
      value <- function() { state$n <- state$n + 1L; bw }
      set.seed(315); before <- .Random.seed
      fit <- switch(shape,
        native=ff(value(),txdat=x,tydat=y,se=FALSE),
        all.positional=ff(value(),x,y,se=FALSE),
        partial=ff(bw=value(),txdat=x,tydat=y,se=FALSE),
        wrapper=(function(...) ff(...))(value(),txdat=x,tydat=y,se=FALSE))
      expect_identical(.Random.seed,before)
      expect_identical(state$n, 1L)
      expect_equal(c(fit$bws$ybw,fit$bws$xbw),bw,tolerance=0)
      expect_equal(fitted(fit),fitted(named),tolerance=0)
    }
  }
})
