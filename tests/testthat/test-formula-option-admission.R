test_that("formula option admission inspects values without executing them", {
  ns <- asNamespace(environmentName(environment(npregbw)))
  admit <- get(".np_formula_options", ns)
  hit <- 0L
  safe <- list(NULL, 1, TRUE, "a", list(quote(hit <<- hit + 1L)),
               expression(hit <<- hit + 1L), function() { hit <<- hit + 1L })
  expect_identical(admit(safe, "test"), safe)
  expect_identical(hit, 0L)
  structural <- list(formula = y ~ x)
  expect_identical(admit(structural, "test"), structural)
  expect_error(admit(list(formula = structure(quote({hit <<- 1L; .4}), class="formula")), "test"), "executable option")
  bad <- list(quote({hit <<- hit + 1L; .4}), as.name("sentinel"),
              compiler::compile(quote({hit <<- hit + 1L; .4})),
              structure(quote({hit <<- hit + 1L; .4}), class = "formula"))
  for (value in bad) {
    expect_error(admit(list(control = value), "test"), "executable option.*control")
    expect_identical(hit, 0L)
  }
  expect_error(admit(list(bad[[1]]), "test"), "#1", fixed = TRUE)
  # Original caller expressions are ordinary evaluated values, including callbacks.
  value <- { hit <- hit + 1L; .4 }
  expect_identical(admit(list(control = value), "test"), list(control = .4))
  expect_identical(hit, 1L)
})

test_that("formula constructors refuse second execution before data can change", {
  old <- options(np.messages = FALSE); on.exit(options(old))
  d <- data.frame(x = seq(-1, 1, length.out = 30), z = cos(seq_len(30)),
                  y = sin(seq_len(30)/4))
  calls <- list(
    reg = list(fun = npregbw, args = list(formula = y~x, data=d, bws=.4)),
    index = list(fun = npindexbw, args = list(formula = y~x+z, data=d, bws=c(1,.5,.4))),
    scoef = list(fun = npscoefbw, args = list(formula = y~x|z, data=d, bws=.4)),
    density = list(fun = npudensbw, args = list(formula = ~x, data=d, bws=.4)),
    distribution = list(fun = npudistbw, args = list(formula = ~x, data=d, bws=.4)),
    cdensity = list(fun = npcdensbw, args = list(formula = y~x, data=d, bws=c(.4,.4))),
    cdistribution = list(fun = npcdistbw, args = list(formula = y~x, data=d, bws=c(.4,.4))))
  for (spec in calls) {
    for (code in list(quote({ydat <- rev(ydat); .1}),
                     quote({mf[[1L]] <- rev(mf[[1L]]); .1}),
                     quote({rm(list="mf"); .1}),
                     quote({assign("tampered", TRUE, envir=capture); .1}))) {
      args <- c(spec$args, list(bandwidth.compute=FALSE,
                    scale.factor.search.lower=code))
      expect_error(do.call(spec$fun, args, quote=TRUE), "executable option")
    }
  }
  # This explicit control is outside ... and must not bypass the same guard.
  expect_error(nplsqregbw(bws=y~x,data=d,tau=quote({mf[[1L]]<-0; .5}),
             bw=.4,bandwidth.compute=FALSE,scale=rep(1,30),delta=1),
             "executable option.*tau")
})

test_that("formula refusal preserves ordinary forcing and retained values", {
  old <- options(np.messages=FALSE);on.exit(options(old))
  d <- data.frame(x=seq(-1,1,length.out=30),y=sin(seq_len(30)/4))
  counts <- c(data=0L, factory=0L, na=0L, control=0L)
  producer <- function(){counts['data']<<-counts['data']+1L;d}
  policy <- function(){counts['factory']<<-counts['factory']+1L;
    function(z){counts['na']<<-counts['na']+1L;stats::na.omit(z)}}
  control <- function(){counts['control']<<-counts['control']+1L;.4}
  set.seed(754);before <- .Random.seed
  a <- npregbw(y~x,data=producer(),na.action=policy(),bws=control(),bandwidth.compute=FALSE)
  expect_identical(unname(counts),rep(1L,4L))
  expect_identical(.Random.seed,before)
  b <- npregbw(y~x,data=d,bws=.4,bandwidth.compute=FALSE)
  expect_identical(a$bw,b$bw)
  expect_equal(fitted(npreg(a)),fitted(npreg(b)),tolerance=0)
  expect_identical(unname(counts),rep(1L,4L))
  # Native/direct data retain their separate successful path.
  native <- npregbw(xdat=d['x'],ydat=d$y,bws=.4,bandwidth.compute=FALSE)
  expect_identical(a$bw,native$bw)
  expect_equal(fitted(npreg(a)),fitted(npreg(native,txdat=d['x'],tydat=d$y)),tolerance=0)
})


test_that("IV formula evaluation values cannot execute at the native handoff", {
  old <- options(np.messages = FALSE); on.exit(options(old))
  d <- data.frame(x = seq(-1, 1, length.out = 30), z = cos(seq_len(30)),
                  y = sin(seq_len(30)/4))
  for (fun in list(npregiv, npregivderiv)) {
    controls <- if (identical(fun, npregiv)) c("zeval", "xeval") else
      c("zeval", "weval", "xeval")
    for (control in controls) {
      args <- list(y = y~x|z, data = d)
      args[[control]] <- quote({stop("SECOND_EXECUTION_OCCURRED"); NULL})
      expect_error(do.call(fun, args, quote = TRUE),
                   paste0("executable option.*", control))
    }
  }
})


test_that("raw formula estimators refuse fit-only executable values before reentry", {
  old <- options(np.messages = FALSE); on.exit(options(old))
  d <- data.frame(x = seq(-1, 1, length.out = 30), z = cos(seq_len(30)),
                  y = sin(seq_len(30)/4), g = factor(rep(c("a", "b"), 15)))
  specs <- list(
    reg = list(fun=npreg, args=alist(y~x, data=d, nmulti=1L), option="gradients"),
    reg.named = list(fun=npreg, args=alist(formula=y~x, data=d, nmulti=1L), option="gradients"),
    reg.fixed = list(fun=npreg, args=alist(y~x, data=d, bws=.4, bandwidth.compute=FALSE), option="gradients"),
    dens = list(fun=npudens, args=alist(~x, data=d, bws=.4, bandwidth.compute=FALSE), option="newdata"),
    dist = list(fun=npudist, args=alist(~x, data=d, bws=.4, bandwidth.compute=FALSE), option="newdata"),
    cdens = list(fun=npcdens, args=alist(y~x, data=d, bws=c(.4,.4), bandwidth.compute=FALSE), option="gradients"),
    cdist = list(fun=npcdist, args=alist(y~x, data=d, bws=c(.4,.4), bandwidth.compute=FALSE), option="gradients"),
    qreg.raw = list(fun=npqreg, args=alist(y~x, data=d, bwmethod="normal-reference"), option="gradients"),
    dist.raw = list(fun=npudist, args=alist(~x, data=d, bwmethod="normal-reference"), option="newdata"),
    qreg = list(fun=npqreg, args=alist(y~x, data=d, bws=c(.4,.4), bandwidth.compute=FALSE), option="gradients"),
    mode = list(fun=npconmode, args=alist(g~x, data=d, bws=c(.1,.4), bandwidth.compute=FALSE), option="proper"),
    scoef = list(fun=npscoef, args=alist(y~x|z, data=d, bws=.4, bandwidth.compute=FALSE), option="residuals"),
    plreg = list(fun=npplreg, args=alist(y~x|z, data=d, bws=matrix(.4,2,1), bandwidth.compute=FALSE), option="residuals"),
    index = list(fun=npindex, args=alist(y~x+z, data=d, bws=c(1,.5,.4), bandwidth.compute=FALSE), option="gradients"))
  hit <- 0L
  for (spec in specs) {
    # An ordinary expression produces an executable value; it is the VALUE
    # that must be refused before do.call can mutate the prepared frame.
    code <- quote(quote({hit <<- hit + 1L; frame.state$frame[[1L]] <- 0; FALSE}))
    args <- spec$args; args[[spec$option]] <- code
    call <- as.call(c(list(spec$fun), args))
    expect_error(eval(call), "executable option", info=spec$option)
    expect_identical(hit, 0L)
  }
})


test_that("every formula transport rejects executable options", {
 ns <- asNamespace(environmentName(environment(npregbw)))
 old <- options(np.messages=FALSE);on.exit(options(old))
 d<-data.frame(x=seq(-1,1,length.out=30),z=cos(seq_len(30)),y=sin(seq_len(30)/4)); dc<-d;dc$y<-factor(rep(c('a','b'),15))
 reg<-npregbw(y~x,data=d,bws=.4,bandwidth.compute=FALSE)
 dens<-npudensbw(~x,data=d,bws=.4,bandwidth.compute=FALSE)
 dist<-npudistbw(~x,data=d,bws=.4,bandwidth.compute=FALSE)
 cd<-npcdensbw(y~x,data=d,bws=c(.4,.4),bandwidth.compute=FALSE)
 cdist<-npcdistbw(y~x,data=d,bws=c(.4,.4),bandwidth.compute=FALSE)
 catbw<-npcdensbw(y~x,data=dc,bws=c(.1,.4),bandwidth.compute=FALSE)
 ix<-npindexbw(y~x+z,data=d,bws=c(1,.5,.4),bandwidth.compute=FALSE)
 sc<-npscoefbw(y~x|z,data=d,bws=.4,bandwidth.compute=FALSE)
 pl<-npplregbw(y~x|z,data=d,bws=matrix(.4,2,1),bandwidth.compute=FALSE)
 state<-new.env(parent=emptyenv());state$frame<-model.frame(y~x,data=d)
 cases<-list(
  npregbw.formula=list(formula=y~x,data=d),npindexbw.formula=list(formula=y~x+z,data=d),
  npscoefbw.formula=list(formula=y~x|z,data=d),npplregbw.formula=list(formula=y~x|z,data=d),
  npudensbw.formula=list(formula=~x,data=d),npudistbw.formula=list(formula=~x,data=d),
  npcdensbw.formula=list(formula=y~x,data=d),npcdistbw.formula=list(formula=y~x,data=d),
  nplsqregbw.formula=list(bws=y~x,data=d),nplsqreg.formula=list(bws=y~x,data=d),
  npreg.formula=list(bws=reg),npindex.formula=list(bws=ix),npscoef.formula=list(bws=sc),
  npplreg.formula=list(bws=pl),npudens.formula=list(bws=dens),npudist.formula=list(bws=dist),
  npcdens.formula=list(bws=cd),npcdist.formula=list(bws=cdist),npconmode.formula=list(bws=catbw),
  npqreg.formula=list(bws=cdist),npreghat.formula=list(bws=reg),npksum.formula=list(formula=y~x,data=d,bws=.4),
  npsigtest.formula=list(bws=reg,.np.formula.state=state),npcopula.formula=list(bws=~x+z,data=d),
  npregiv.formula=list(y=y~x|z,data=d),npregivderiv.formula=list(y=y~x|z,data=d))
 z<-list();set.seed(754)
 for(owner in names(cases)){
  args<-cases[[owner]];arg<-if(owner=='npksum.formula')'kernel.pow' else if(owner=='npqreg.formula')'tol' else 'r75probe'
  code<-quote({stop('SECOND_EXECUTION_OCCURRED');1});args[[arg]]<-code
  before<-.Random.seed;err<-tryCatch({do.call(get(owner,ns),args,quote=TRUE);'NO_ERROR'},error=conditionMessage)
  ok<-grepl('executable option',err,fixed=TRUE);z[[owner]]<-data.frame(owner,option=arg,refused=ok,rng_same=identical(.Random.seed,before),message=err)
  expect_true(ok, info=paste(owner,err)); expect_identical(.Random.seed,before)
 }

})


test_that("explicit formula reentry controls use the same refusal", {
  old <- options(np.messages=FALSE);on.exit(options(old))
  d <- data.frame(x=seq(-1,1,length.out=30),z=cos(seq_len(30)),
                  y=sin(seq_len(30)/4),g=factor(rep(c("a","b"),15)))
  for (control in c("level","proper","proper.control")) {
    args <- alist(g~x,data=d,bws=c(.1,.4),bandwidth.compute=FALSE)
    args[[control]] <- quote(quote({stop("SECOND_EXECUTION_OCCURRED");NULL}))
    expect_error(eval(as.call(c(list(npconmode),args))), "executable option")
  }
  for (control in c("u","n.quasi.inv","er.quasi.inv")) {
    args <- alist(~x+z,data=d,bws=c(.4,.4),bandwidth.compute=FALSE,neval=3)
    args[[control]] <- quote(quote({stop("SECOND_EXECUTION_OCCURRED");NULL}))
    expect_error(eval(as.call(c(list(npcopula),args))), "executable option")
  }
  # The special single-index formula plus explicit native-data reentry also
  # prepares a frame before the default method's constructor and fit calls.
  expect_error(npindex(bws=y~x+z,txdat=d[c("x","z")],tydat=d$y,
    nmulti=1L,
    gradients=quote({stop("SECOND_EXECUTION_OCCURRED");FALSE})),"executable option")
})


test_that("raw formula significance-test options cannot execute in the prepared owner", {
  old <- options(np.messages=FALSE);on.exit(options(old))
  d <- data.frame(x=seq(-1,1,length.out=30),y=sin(seq_len(30)/4))
  expect_error(npsigtest(bws=y~x,data=d,nmulti=1L,B=9L,
    pivot=quote({stop("SECOND_EXECUTION_OCCURRED");FALSE})),
    "executable option.*pivot")
})
