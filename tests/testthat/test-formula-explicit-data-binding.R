test_that("least-squares quantile formula owners bind wrapper data values", {
  old <- options(np.messages=FALSE)
  on.exit(options(old),add=TRUE)
  fixture <- data.frame(x=seq(.05,.95,length.out=24),
    y=sin(seq_len(24)/3)+seq_len(24)/24)
  fo <- y~x
  dd <- transform(fixture,y=y+100)
  controls <- list(scale=rep(1,24),bandwidth.compute=FALSE,nomad=FALSE,
    regtype="ll",nmulti=1L,itmax=5L,random.seed=17L)
  for(fun in list(nplsqregbw,nplsqreg)) {
    reference <- do.call(fun,c(list(bws=fo,data=fixture),controls))
    wrapped <- function(dd) fun(fo,data=dd,scale=rep(1,nrow(dd)),
      bandwidth.compute=FALSE,nomad=FALSE,regtype="ll",nmulti=1L,
      itmax=5L,random.seed=17L)
    expression <- function(dd) fun(fo,data=dd[,],scale=rep(1,nrow(dd)),
      bandwidth.compute=FALSE,nomad=FALSE,regtype="ll",nmulti=1L,
      itmax=5L,random.seed=17L)
    for(actual in list(wrapped(fixture),expression(fixture))) {
      expect_identical(class(actual),class(reference))
      if(inherits(actual,"lsqregression")) {
        expect_identical(fitted(actual),fitted(reference))
      } else {
        expect_identical(actual$xdat,reference$xdat)
        expect_identical(actual$ydat,reference$ydat)
        expect_identical(actual$reg.bws$bw,reference$reg.bws$bw)
      }
    }
  }
})

test_that("formula owners bind explicit data values without rebinding their names", {

  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  fixture <- data.frame(x = seq(.1, .9, length.out = 24),
                        z = sin(1:24), y = cos(1:24))
  fixture$cl <- factor(rep(c("a", "b"), 12))
  dd <- transform(fixture, x = rev(x), y = y + 10)
  build <- function(fun, formula, ...) do.call(fun, c(list(
    formula = formula, data = fixture, bandwidth.compute = FALSE), list(...)))
  rb <- build(npregbw, y ~ x, bws = .25, regtype = "ll")
  db <- build(npudensbw, ~ x, bws = .25)
  ub <- build(npudistbw, ~ x, bws = .25)
  cb <- build(npcdensbw, y ~ x, bws = c(.3, .25))
  qb <- build(npcdistbw, y ~ x, bws = c(.3, .25))
  mb <- build(npcdensbw, cl ~ x, bws = c(.2, .25))
  sb <- build(npscoefbw, y ~ x | z, bws = .3, regtype = "lc")
  pb <- build(npplregbw, y ~ x | z, bws = matrix(.3, 2, 1), regtype = "lc")
  ib <- build(npindexbw, y ~ x + z, bws = c(1, .2, .3), regtype = "lc")
  cases <- list(npreg = rb, npudens = db, npudist = ub, npcdens = cb,
    npcdist = qb, npqreg = qb, npconmode = mb, npscoef = sb, npplreg = pb,
    npindex = ib, npreghat = rb, npsigtest = rb)
  fields <- c("mean", "merr", "grad", "gerr", "dens", "derr", "dist", "disterr",
    "condens", "conderr", "condist", "quantile", "quanterr", "conmode",
    "probabilities", "beta", "In", "P")
  payload <- function(obj) {
    if (is.matrix(obj) || is.atomic(obj)) {
      attr(obj, "call") <- NULL
      return(obj)
    }
    ans <- obj[intersect(names(obj), fields)]
    expect_gt(length(ans), 0L)
    ans
  }
  for (name in names(cases)) {
    fun <- get(name, mode = "function")
    bw <- cases[[name]]
    extra <- if (name == "npsigtest") list(B = 9, random.seed = 91) else list()
    reference <- do.call(fun, c(list(bws = bw, data = fixture), extra))
    local.call <- function(dd) {
      if (name == "npsigtest") return(fun(bw, data = dd, B = 9, random.seed = 91))
      fun(bw, data = dd)
    }
    expression.call <- function(dd) {
      if (name == "npsigtest") return(fun(bw, data = dd[, , drop = FALSE], B = 9, random.seed = 91))
      fun(bw, data = dd[, , drop = FALSE])
    }
    expect_identical(payload(local.call(fixture)), payload(reference), info = name)
    expect_identical(payload(expression.call(fixture)), payload(reference), info = name)
    ## Retained-call behavior is untouched when no override is supplied.
    expect_identical(payload(do.call(fun, c(list(bws = bw), extra))),
                     payload(reference), info = name)
  }
})

test_that("explicit regression data preserves formula scope, subset, NA and prediction", {

  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  formula.env <- new.env(parent = environment())
  formula.env$offset.value <- .5
  formula.env$shift <- eval(quote(function(x) x + offset.value), formula.env)
  fo <- as.formula("y ~ shift(x)", env = formula.env)
  fixture <- data.frame(y = sin(1:24), x = (1:24) / 24)
  fixture$y[4] <- NA_real_
  fixture$keep <- seq_len(nrow(fixture)) > 2
  bw <- npregbw(fo, data = fixture, subset = keep, na.action = na.omit,
                bws = .3, bandwidth.compute = FALSE, regtype = "ll")
  reference <- do.call(npreg, list(bws = bw, data = fixture, se = TRUE, gradients = TRUE))
  local.call <- function(dd) {
    offset.value <- 900
    npreg(bw, data = dd, se = TRUE, gradients = TRUE)
  }
  actual <- local.call(fixture)
  for (field in c("mean", "merr", "grad", "gerr", "rows.omit", "nobs.omit"))
    expect_identical(actual[[field]], reference[[field]])
  nd <- data.frame(x = c(.2, .6))
  expect_identical(predict(actual, newdata = nd), predict(reference, newdata = nd))
  expect_identical(se(actual), se(reference))
  expect_identical(gradients(actual), gradients(reference))
  expect_output(summary(actual))
  expect_error(predict(actual, newdata = data.frame(wrong = 1:2)), "newdata")
})
