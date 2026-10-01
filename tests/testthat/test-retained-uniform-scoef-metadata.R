# This body is copied byte-for-byte into each package's new test file.
test_that("smooth-coefficient reconstruction distinguishes retained uniform order", {
  old <- options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  pkg <- getNamespaceName(environment(npscoefbw))
  ns <- asNamespace(pkg)
  context <- get(".np_retained_uniform_constructor",ns)
  collect <- function(expr) {
    warnings <- character()
    value <- withCallingHandlers(expr,warning=function(w) {
      warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")
    })
    list(value=value,warnings=warnings)
  }
  message <- unname(get(".np_io_prefix_text",ns)(
    "ignoring kernel order specified with uniform kernel type"))
  constructor <- function(ckertype, order, payload) {
    warning(message,call.=FALSE)
    list(ckertype=ckertype,order=order,payload=payload)
  }
  payload <- list(x=c(-1,0,1),order=4L)
  a <- collect(context(constructor,ckertype="uniform",order=4L,payload=payload))
  expect_identical(a$value,list(ckertype="uniform",order=4L,payload=payload))
  expect_length(a$warnings,0L)
  expect_identical(collect(context(constructor,ckertype="gaussian",order=4L,
    payload=payload))$warnings,message)
  expect_identical(collect(context(function(...)warning("unrelated",call.=FALSE),
    ckertype="uniform"))$warnings,"unrelated")
  expect_error(context(function(...)stop("unrelated-error"),ckertype="uniform"),
    "unrelated-error",fixed=TRUE)

  n <- 60L
  x <- data.frame(x=seq(-1,1,length.out=n))
  z <- data.frame(z=sin(seq_len(n)))
  y <- x$x*(1+z$z)+cos(seq_len(n))/10
  for(kernel in c("uniform","gaussian"))for(explicit in c(FALSE,TRUE)) {
    args <- list(xdat=x,zdat=z,ydat=y,bws=.6,ckertype=kernel,
      bandwidth.compute=FALSE)
    if(explicit)args$ckerorder <- 4L
    got <- collect(do.call(npscoefbw,args))
    expect_length(got$warnings,as.integer(kernel=="uniform"&&explicit))
    if(length(got$warnings))expect_identical(got$warnings,message)
    expect_identical(as.integer(got$value$ckerorder),if(explicit)4L else 2L)
    expect_equal(as.double(got$value$bw),.6,tolerance=0)
  }
  got <- collect(npscoefbw(xdat=x,zdat=z,ydat=y,bws=.6,ckertype="uniform",
    nmulti=1L,optim.maxit=3L,random.seed=42L))
  expect_length(got$warnings,0L)
  expect_true(is.finite(got$value$fval))
})
