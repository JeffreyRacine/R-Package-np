test_that("CDF single-evaluation admission is opt-in and preserves default results", {
  old <- options(np.messages=FALSE,np.tree=FALSE,np.extendednn=FALSE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(x1=c(rep(0,32),sin((1:16)*sqrt(2))+(1:16)/32),
    x2=c(rep(0,32),sin((1:16)*sqrt(3))+(1:16)/48))
  y <- data.frame(y=sin((1:48)/3)+(1:48)/90)
  fields <- c("objective","num.feval","num.feval.fast")
  for(type in c("fixed","generalized_nn","adaptive_nn"))
    for(degree in list(c(0L,0L),c(1L,1L),c(1L,0L)))
      for(k in if(type=="fixed")1 else c(2,if(type=="adaptive_nn")46 else 47)) {
        bw <- npcdistbw(xdat=x,ydat=y,bwtype=type,regtype="lp",degree=degree,
          nomad=FALSE,bws=rep(k,3L),bandwidth.compute=FALSE)
        args <- list(xdat=x,ydat=y,bws=bw,ngrid=7L,invalid.penalty="baseline")
        default <- do.call(.npcdistbw_eval_only,args)
        explicit.false <- do.call(.npcdistbw_eval_only,c(args,list(return.admissible=FALSE)))
        status <- do.call(.npcdistbw_eval_only,c(args,list(return.admissible=TRUE)))
        args$invalid.penalty <- "dbmax"
        raw <- do.call(.npcdistbw_eval_only,args)
        expect_identical(names(default),fields)
        expect_identical(explicit.false,default)
        expect_identical(status[fields],default)
        expect_identical(status$admissible,.np_nn_raw_objective_valid(raw$objective))
        expect_identical(status$admissible,!(type!="fixed"&&k==2))
        expect_identical(default$num.feval,1L)
      }
})
