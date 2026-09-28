test_that("retained uniform hats preserve quiet metadata and fresh request advisories", {
  old <- options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  pkg <- getNamespaceName(environment(npregbw))
  collect <- function(expr) {
    messages <- classes <- list()
    value <- withCallingHandlers(expr,warning=function(w) {
      messages[[length(messages)+1L]] <<- conditionMessage(w)
      classes[[length(classes)+1L]] <<- class(w)
      invokeRestart("muffleWarning")
    })
    list(value=value,messages=unlist(messages),classes=classes)
  }
  x <- data.frame(x=seq(-1,1,length.out=37L))
  y <- sin(2*x$x)+seq_len(nrow(x))/100
  ex <- x[c(3L,11L,24L),,drop=FALSE]
  for(type in c("fixed","generalized_nn","adaptive_nn")) {
    b <- npregbw(xdat=x,ydat=y,bws=if(type=="fixed") .7 else 12,
      bwtype=type,regtype="lc",ckertype="uniform",bandwidth.compute=FALSE)
    for(eval in list(NULL,ex)) {
      args <- list(bws=b,txdat=x)
      if (!is.null(eval)) args$exdat <- eval
      fit <- collect(do.call(npreg,c(args,list(tydat=y,gradients=TRUE))))
      expect_null(fit$messages)
      for(s in 0:1) {
        h <- collect(do.call(npreghat,c(args,list(s=s,output="matrix"))))
        expect_null(h$messages)
        expect_true(all(is.finite(h$value)))
      }
    }
  }
  yf <- data.frame(y=y)
  for(cdf in c(FALSE,TRUE)) for(type in c("fixed","adaptive_nn")) {
    make <- if(cdf) npcdistbw else npcdensbw
    hat <- if(cdf) npcdisthat else npcdenshat
    b <- make(xdat=x,ydat=yf,bws=rep(if(type=="fixed") .7 else 12,2L),
      bwtype=type,regtype="lc",cxkertype="uniform",cykertype="uniform",bandwidth.compute=FALSE)
    for(s in 0:1) {
      h <- collect(hat(b,txdat=x,tydat=yf,exdat=ex,eydat=yf[c(3L,11L,24L),,drop=FALSE],s=s,output="matrix"))
      a <- collect(hat(b,txdat=x,tydat=yf,exdat=ex,eydat=yf[c(3L,11L,24L),,drop=FALSE],s=s,y=seq_len(37),output="apply"))
      expect_null(h$messages); expect_null(a$messages)
      expect_equal(drop(h$value %*% seq_len(37)),as.vector(a$value),tolerance=2e-11)
    }
    for(role in c("x","y")) {
      fun <- getFromNamespace(if(role=="x") ".npcdhat_make_xkbw" else ".npcdhat_make_ybw",pkg)
      kb <- collect(fun(b,if(role=="x") x else yf))
      expect_null(kb$messages)
      expect_identical(kb$value$ckerorder,b[[paste0("c",role,"kerorder")]])
    }
  }
  for(order in c(2L,4L)) {
    z <- collect(npregbw(xdat=x,ydat=y,bws=.7,regtype="lc",ckertype="uniform",ckerorder=order,bandwidth.compute=FALSE))
    expect_length(z$messages,1L)
    expect_match(z$messages,"ignoring kernel order specified with uniform kernel type",fixed=TRUE)
    expect_true("warning" %in% z$classes[[1L]])
  }
  make <- getFromNamespace(".npcdhat_make_xkbw",pkg)
  bad <- b; bad$cxkerbound <- "invalid"
  expect_error(make(bad,x),"arg")
  expect_null(collect(make(b,x))$messages)
})
