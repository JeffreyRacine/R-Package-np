test_that("conditional plot builders keep retained uniform metadata quiet", {
  old <- options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  pkg <- getNamespaceName(environment(npcdensbw))
  ns <- asNamespace(pkg)
  collect <- function(expr) {
    warnings <- character()
    value <- withCallingHandlers(expr,warning=function(w) {
      warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")
    })
    list(value=value,warnings=warnings)
  }
  x <- data.frame(x=seq(-1,1,length.out=31))
  y <- data.frame(y=sin(x$x)+seq_len(31)/100)
  for(cdf in c(FALSE,TRUE)) for(type in c("fixed","generalized_nn","adaptive_nn"))
    for(kernel in c("uniform","gaussian")) {
      make <- if(cdf)npcdistbw else npcdensbw
      b <- make(xdat=x,ydat=y,bws=rep(if(type=="fixed").7 else 12,2),
        bwtype=type,regtype="lc",cxkertype=kernel,cykertype=kernel,
        bandwidth.compute=FALSE)
      for(role in c("x","xy")) {
        builder <- get(paste0(".np_con_make_kbandwidth_",role),ns)
        args <- list(bws=b,xdat=x);if(role=="xy")args$ydat<-y
        got <- collect(do.call(builder,args))
        expect_length(got$warnings,0L)
        expect_identical(got$value$ckerorder,b$cxkerorder)
        # Test-only lexical substitution restores exactly the original
        # constructor while leaving every argument and numerical value intact.
        original <- builder;env <- new.env(parent=ns)
        env$.npcdhat_retained_kbandwidth <- get("kbandwidth.numeric",ns)
        environment(original) <- env
        before <- collect(do.call(original,args))
        expect_identical(got$value,before$value)
        expect_length(before$warnings,as.integer(kernel=="uniform"))
      }
    }
  b <- npcdensbw(xdat=x,ydat=y,bws=c(.7,.7),regtype="lc",
    cxkertype="uniform",cykertype="uniform",bandwidth.compute=FALSE)
  b$cxkerbound <- "invalid"
  expect_error(get(".np_con_make_kbandwidth_x",ns)(b,x),"arg")
  fresh <- collect(npcdensbw(xdat=x,ydat=y,bws=c(.7,.7),regtype="lc",
    cxkertype="uniform",cxkerorder=4L,bandwidth.compute=FALSE))
  expect_length(fresh$warnings,1L)
  expect_match(fresh$warnings,"ignoring cxkerorder",fixed=TRUE)
})
