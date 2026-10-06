diagnostics_contract_fixture <- function() {
  set.seed(106L)
  d <- data.frame(x = rnorm(36), z = rnorm(36))
  d$y <- as.double(d$x + .2*d$z + rnorm(36) > 0)
  d
}

diagnostics_contract_pool <- function() {
  if ("npRmpi" %in% loadedNamespaces() && !spawn_mpi_slaves())
    skip("MPI pool unavailable in this test context")
}

test_that("single-index diagnostics select the requested sample without changing predictions", {
  diagnostics_contract_pool()
  if ("npRmpi" %in% loadedNamespaces())
    on.exit(close_mpi_slaves(force=TRUE), add=TRUE)
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  d <- diagnostics_contract_fixture(); e <- d[1:7, ]; e$y <- 1-e$y
  bw <- npindexbw(y ~ x+z, data=d, method="kleinspady",
                   bws=c(1,.2,.7), bandwidth.compute=FALSE)
  fit <- npindex(bw, se=FALSE)
  scored <- npindex(bw, newdata=e, se=FALSE)
  unscored <- npindex(bw, newdata=e[c("x","z")], se=FALSE)
  explicit <- npindex(bw, newdata=e, y.eval=FALSE, se=FALSE)
  expect_identical(fitted(scored), fitted(unscored))
  expect_identical(fitted(scored), fitted(explicit))
  expect_identical(scored$diagnostics.sample, "evaluation")
  expect_identical(scored$diagnostics.nobs, 7L)
  expect_equal(sum(scored$confusion.matrix), 7)
  expect_equal(scored$CCR.overall, mean(e$y == pmin(1L,pmax(0L,round(fitted(scored))))))
  expect_identical(unscored$confusion.matrix, NA)
  expect_identical(unscored$diagnostics.sample, "unavailable")
  expect_identical(explicit$confusion.matrix, fit$confusion.matrix)
  expect_output(summary(explicit), "Training diagnostics", fixed=TRUE)
  expect_output(summary(scored), "Evaluation diagnostics", fixed=TRUE)
  expect_output(summary(unscored), "Evaluation diagnostics unavailable", fixed=TRUE)
  e$y[c(2,5)] <- NA_real_
  partial <- npindex(bw, newdata=e, se=FALSE)
  expect_identical(fitted(partial), fitted(scored))
  expect_identical(partial$diagnostics.nobs, 5L)
  expect_equal(sum(partial$confusion.matrix), 5)
  expect_error(npindex(bw,newdata=e[c("x","z")],y.eval=TRUE,se=FALSE), "columns|requires evaluation outcomes")
  expect_error(npindex(bw,newdata=e,y.eval=NA,se=FALSE), "y.eval", fixed=TRUE)
})

test_that("continuous index diagnostics use transformed outcomes only for scoring", {
  diagnostics_contract_pool()
  if ("npRmpi" %in% loadedNamespaces())
    on.exit(close_mpi_slaves(force=TRUE), add=TRUE)
  old <- options(np.messages=FALSE); on.exit(options(old),add=TRUE)
  d <- diagnostics_contract_fixture(); d$y <- exp(d$x + .2*d$z)
  b <- npindexbw(log(y) ~ x+z, data=d, bws=c(1,.2,.7), bandwidth.compute=FALSE)
  e <- d[1:7, ]; e$y <- e$y*2
  f <- npindex(b,newdata=e,se=FALSE)
  x <- npindex(b,newdata=e[c("x","z")],se=FALSE)
  expect_identical(fitted(f),fitted(x))
  expect_equal(f$MSE,mean((log(e$y)-fitted(f))^2))
  expect_true(is.na(x$MSE))
  e$y[] <- NA_real_
  missing <- npindex(b,newdata=e,se=FALSE)
  expect_identical(fitted(missing),fitted(x))
  expect_identical(missing$diagnostics.nobs,0L)
  expect_true(is.na(missing$MSE))
})

test_that("conditional mode formula and native newdata share outcome diagnostics", {
  diagnostics_contract_pool()
  if ("npRmpi" %in% loadedNamespaces())
    on.exit(close_mpi_slaves(force=TRUE), add=TRUE)
  old <- options(np.messages=FALSE); on.exit(options(old),add=TRUE)
  d <- data.frame(x=seq(-1,1,length.out=36),
                  y=factor(rep(c("2","5","9"),12),levels=c("9","2","5")))
  e <- d[c(2,8,11,16,20,25,31), ]; e$y <- factor(rep("5",7),levels=levels(d$y))
  for (formula in c(TRUE,FALSE)) {
    b <- if (formula) npcdensbw(y~x,data=d,bws=c(.2,.5),bandwidth.compute=FALSE) else
      npcdensbw(xdat=d["x"],ydat=d["y"],bws=c(.2,.5),bandwidth.compute=FALSE)
    s <- npconmode(b,newdata=e,probabilities=TRUE)
    x <- npconmode(b,newdata=e["x"],probabilities=TRUE)
    expect_identical(s$conmode,x$conmode)
    expect_identical(s$probabilities,x$probabilities)
    expect_equal(sum(s$confusion.matrix),7)
    expect_equal(s$CCR.overall,mean(as.character(s$conmode)==as.character(e$y)))
    expect_identical(x$confusion.matrix,NA)
    expect_identical(s$diagnostics.sample,"evaluation")
    expect_output(summary(x),"Evaluation diagnostics unavailable",fixed=TRUE)
    partial <- e; partial$y[c(1,3)] <- NA
    p <- npconmode(b,newdata=partial,probabilities=TRUE)
    expect_identical(p$probabilities,s$probabilities)
    expect_identical(p$diagnostics.nobs,5L)
    expect_equal(sum(p$confusion.matrix),5)
  }
})

test_that("evaluation row maps and native precedence do not change fitted quantities", {
  diagnostics_contract_pool()
  if ("npRmpi" %in% loadedNamespaces())
    on.exit(close_mpi_slaves(force=TRUE), add=TRUE)
  old <- options(np.messages=FALSE); on.exit(options(old),add=TRUE)
  set.seed(10607)
  d <- data.frame(x=rnorm(80),z=rnorm(80))
  d$y <- as.double(d$x+.2*d$z+rnorm(80)>0)
  y <- d$y # Native metadata retains this unambiguous response name.
  e <- d[1:9, ]; e$y <- 1-e$y; e$x[3] <- NA_real_; e$y[5] <- NA_real_
  for (formula in c(TRUE,FALSE)) {
    b <- if(formula) npindexbw(y~x+z,data=d,method="kleinspady",
        bws=c(1,.2,.7),bandwidth.compute=FALSE,na.action=na.exclude) else
      npindexbw(xdat=d[c("x","z")],ydat=y,method="kleinspady",
        bws=c(1,.2,.7),bandwidth.compute=FALSE)
    eval.data <- e
    if (!formula) eval.data[[b$ynames]] <- e$y
    rng <- .Random.seed
    a <- npindex(b,newdata=eval.data,se=TRUE,gradients=TRUE)
    x <- npindex(b,newdata=e[c("x","z")],se=TRUE,gradients=TRUE)
    if (!("npRmpi" %in% loadedNamespaces())) expect_identical(.Random.seed,rng)
    expect_identical(a$bws$beta,b$beta)
    expect_identical(fitted(a),fitted(x))
    expect_identical(se(a),se(x))
    expect_identical(gradients(a),gradients(x))
    expect_identical(a$diagnostics.nobs,7L)
    expect_equal(sum(a$confusion.matrix),7)
    explicit <- npindex(b,newdata=e,eydat=e$y,se=TRUE,gradients=TRUE,y.eval=FALSE)
    expect_identical(explicit$confusion.matrix,a$confusion.matrix)
    native <- npindex(b,exdat=e[c("x","z")],eydat=e$y,se=TRUE,gradients=TRUE)
    expect_identical(fitted(native),fitted(a))
    expect_identical(native$confusion.matrix,a$confusion.matrix)
    expect_equal(predict(a,newdata=e[c("x","z")],se.fit=TRUE,gradients=TRUE)$fit,fitted(x),tolerance=0)
  }
  d$y <- factor(ifelse(d$y==1,"5","2"),levels=c("5","2"))
  e <- d[1:9, ]; e$x[3] <- NA_real_; e$y[5] <- NA
  for (formula in c(TRUE,FALSE)) {
    b <- if(formula) npcdensbw(y~x+z,data=d,bws=c(.2,.7,.7),
        bandwidth.compute=FALSE,na.action=na.exclude) else
      npcdensbw(xdat=d[c("x","z")],ydat=d["y"],bws=c(.2,.7,.7),bandwidth.compute=FALSE)
    a <- npconmode(b,newdata=e,probabilities=TRUE,se=TRUE,gradients=TRUE)
    x <- npconmode(b,newdata=e[c("x","z")],probabilities=TRUE,se=TRUE,gradients=TRUE)
    expect_identical(a$probabilities,x$probabilities)
    expect_identical(a$probability.errors,x$probability.errors)
    expect_identical(a$probability.gradients,x$probability.gradients)
    expect_equal(sum(a$confusion.matrix),7)
    explicit <- npconmode(b,newdata=e,eydat=e["y"],probabilities=TRUE)
    expect_identical(explicit$confusion.matrix,a$confusion.matrix)
    native <- npconmode(b,exdat=e[c("x","z")],eydat=e["y"],probabilities=TRUE)
    expect_identical(native$conmode,a$conmode)
    expect_identical(native$confusion.matrix,a$confusion.matrix)
    expect_identical(predict(a,newdata=e[c("x","z")],type="prob",se.fit=TRUE,gradients=TRUE)$fit,x$probabilities)
  }
})
