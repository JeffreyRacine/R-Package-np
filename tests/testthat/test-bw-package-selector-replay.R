test_that("MPI selector replay retains its package owner under global masking", {
  skip_on_cran()
  lines <- c(
    "library(npRmpi)",
    "npRmpi.init(nslaves=1L,quiet=TRUE)",
    "on.exit(npRmpi.quit(force=TRUE),add=TRUE)",
    "options(np.messages=FALSE)",
    "ns<-asNamespace('npRmpi')",
    "selectors<-c('npudensbw','npudistbw','npcdensbw','npcdistbw','npregbw','npscoefbw','npindexbw','npplregbw')",
    "assign('.selector_owner_hits',setNames(integer(length(selectors)),selectors),.GlobalEnv)",
    "for(n in selectors)trace(n,where=ns,print=FALSE,tracer=bquote(.GlobalEnv$.selector_owner_hits[[.(n)]] <- .GlobalEnv$.selector_owner_hits[[.(n)]]+1L))",
    "on.exit(for(n in selectors)untrace(n,where=ns),add=TRUE)",
    "assign('.selector_mask_hits',setNames(integer(length(selectors)),selectors),.GlobalEnv)",
    "for(n in selectors)assign(n,local({name<-n;function(...){.GlobalEnv$.selector_mask_hits[[name]] <- .GlobalEnv$.selector_mask_hits[[name]]+1L;get(name,ns)(...)}}),.GlobalEnv)",
    "set.seed(314);x<-data.frame(x=runif(32));z<-data.frame(z=runif(32));y<-sin(3*x$x)+cos(1:32)*.15;dy<-data.frame(y=y)",
    "calls<-list(npudensbw=quote(npRmpi::npudens(tdat=x,nmulti=1L,itmax=20L)), npudistbw=quote(npRmpi::npudist(tdat=x,nmulti=1L,itmax=20L)), npcdensbw=quote(npRmpi::npcdens(txdat=x,tydat=dy,nmulti=1L,itmax=20L)), npcdistbw=quote(npRmpi::npcdist(txdat=x,tydat=dy,nmulti=1L,itmax=20L)), npregbw=quote(npRmpi::npreg(txdat=x,tydat=y,nmulti=1L,itmax=20L)), npscoefbw=quote(npRmpi::npscoef(txdat=x,tydat=y,tzdat=z,nmulti=1L,optim.maxit=20L)), npindexbw=quote(npRmpi::npindex(txdat=x,tydat=y,nmulti=1L,itmax=20L)), npplregbw=quote(npRmpi::npplreg(txdat=x,tydat=y,tzdat=z,nmulti=1L,itmax=20L)))",
    "for(n in names(calls)){before<-.GlobalEnv$.selector_owner_hits[[n]];fit<-eval(calls[[n]]);v<-fitted(fit);if(.GlobalEnv$.selector_mask_hits[[n]]!=0L)stop('global selector mask was evaluated');stopifnot(length(v)>0L,all(is.finite(v)),.GlobalEnv$.selector_owner_hits[[n]]>before);cat('SELECTOR_REPLAY_OWNER',n,'\\n')}",
    "cat('SELECTOR_REPLAY_GLOBAL_MASK_PASS\\n')"
  )
  result <- npRmpi_run_isolated_contract(lines,
    marker="SELECTOR_REPLAY_GLOBAL_MASK_PASS",timeout=120L)
  if(is.null(result))skip("installed subprocess library unavailable")
  expect_identical(result$status,0L,info=paste(result$output,collapse="\n"))
  expect_true(result$witnessed)
})
