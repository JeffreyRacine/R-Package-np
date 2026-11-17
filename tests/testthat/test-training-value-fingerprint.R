h_function <- function(name) {
  get(name, environment(getFromNamespace(".npRmpi_autodispatch_fingerprint", "npRmpi")),
      inherits = TRUE)
}

test_that("training certificates describe values independently of compact storage", {
  cert <- h_function(".npRmpi_training_value_certificate")
  plain <- function(x) {
    a <- attributes(x)
    v <- if (is.object(x)) unclass(x) else x
    out <- if (is.list(v)) lapply(v, plain) else .subset(v, seq_along(v))
    attributes(out) <- lapply(a, plain)
    out
  }
  cases <- list(1:40000, as.double(1:40000), as.character(1:40000),
    sort(sin(1:40000)), c(TRUE, FALSE, NA), c(NA_real_, NaN, Inf, -Inf, -0, 0),
    c(1+2i, NA_complex_), as.raw(0:255), matrix(1:12, ncol=3),
    factor(c("b", "a", NA)), ordered(c("b", "a")),
    as.Date(1:10, origin="2000-01-01"), as.POSIXct(1:10, origin="1970-01-01"),
    as.difftime(1:10, units="secs"), data.frame(x=1:40000, y=as.character(1:40000)))
  for (x in cases) {
    expect_false(is.null(cert(x)))
    expect_identical(cert(x), cert(plain(x)))
    expect_identical(cert(x), cert(unserialize(serialize(x, NULL))))
  }
  expect_null(cert(NULL))
  expect_null(cert(structure(1:3, class="custom")))
  expect_null(cert(list(callback=function(x) x)))
  expect_null(cert(new.env(parent=emptyenv())))
})

test_that("training certificates retain type, attributes and all data mutations", {
  cert <- h_function(".npRmpi_training_value_certificate")
  expect_false(identical(cert(1L), cert(1)))
  expect_false(identical(cert(0), cert(-0)))
  expect_false(identical(cert(NA_real_), cert(NaN)))
  expect_false(identical(cert(structure(1, a=1, b=2)), cert(structure(1, b=2, a=1))))
  d <- data.frame(x=1:40000, y=factor(rep(c("a", "b"),20000)))
  reference <- cert(d)
  for (i in c(1L, 32768L, 32769L, 40000L)) {
    changed <- d; changed$x[i] <- -1L
    expect_false(identical(reference, cert(changed)))
  }
  changed <- d; row.names(changed)[1L] <- "new"
  expect_false(identical(reference, cert(changed)))
  changed <- d; levels(changed$y) <- c("c", "b")
  expect_false(identical(reference, cert(changed)))
  changed <- d; attr(changed$x, "units") <- "new"
  expect_false(identical(reference, cert(changed)))
})

test_that("projection is restricted to named training roles and exact containers", {
  project <- h_function(".npRmpi_training_value_projection")
  d <- data.frame(x=1:10)
  x <- structure(list(call=quote(npreg(xdat=d, ydat=, zdat=NULL)),
    .np.native.training=list(xdat=d, ydat=NULL, zdat=NULL),
    .np.formula.training=list(frame=NULL, subset=quote(x>0)),
    xdati=list(n=10L)), class="rbandwidth")
  p <- project(x)
  expect_true(p$changed)
  expect_identical(p$value$call, x$call)
  expect_identical(p$value$.np.native.training$ydat, NULL)
  expect_identical(p$value$.np.formula.training, x$.np.formula.training)
  expect_identical(p$value$xdati, x$xdati)
  expect_identical(x$.np.native.training$xdat, d)
  for (cl in c("rbandwidth", "lsqregressionbandwidth", "plbandwidth", "sibandwidth",
               "scbandwidth", "conbandwidth", "condbandwidth", "bandwidth", "dbandwidth")) {
    class(x) <- cl; expect_true(project(x)$changed)
  }
  for (cl in c("npregression", "copula", "user.bandwidth")) {
    class(x) <- cl; expect_identical(project(x), list(value=x, changed=FALSE))
  }
  class(x) <- c("user", "rbandwidth")
  expect_identical(project(x), list(value=x, changed=FALSE))
  expect_false(project(pairlist(.np.native.training=list(xdat=d)))$changed)
  bare <- substitute(.np.literal(value=VALUE), list(VALUE=d))
  expect_identical(project(bare), list(value=bare, changed=FALSE))
  call <- substitute(npreg(xdat=VALUE, ydat=, zdat=NULL),list(VALUE=d))
  p <- project(list(call=call))
  expect_true(p$changed)
  expect_identical(names(p$value$call), names(call))
  expect_identical(p$value$call[["zdat"]], NULL)
  expect_identical(p$value$call[["ydat"]], quote(expr=))
})

test_that("new value reads bypass class methods", {
  cert <- h_function(".npRmpi_training_value_certificate")
  # Private methods in the helper's lexical chain catch accidental generic reads.
  env <- new.env(parent=environment(cert))
  count <- 0L
  trap <- function(...) { count <<- count+1L; stop("unexpected class dispatch") }
  for (cl in c("Date", "factor", "data.frame"))
    for (generic in c("length", "[", "[[", "as.list"))
      assign(paste0(generic,".",cl), trap, env)
  f <- cert; environment(f) <- env
  values <- list(as.Date(1:3, origin="2000-01-01"),factor(c("a","b")),data.frame(x=1:3))
  for (x in values) expect_identical(f(x), cert(x))
  expect_identical(count, 0L)
})

test_that("projection failure preserves the incumbent uncertifiable result", {
  fp <- h_function(".npRmpi_autodispatch_fingerprint")
  env <- new.env(parent=environment(fp))
  env$.npRmpi_training_value_projection <- function(x) stop("forced certificate failure")
  probe <- fp; environment(probe) <- env
  x <- list(.np.native.training=list(xdat=data.frame(x=1:3)))
  expect_identical(probe(x), NA_character_)
  expect_identical(probe(list(a=1), policy="context.v1"),
                   fp(list(a=1), policy="context.v1"))
})

test_that("unprojected objects retain the incumbent complete fingerprint", {
  fp <- h_function(".npRmpi_autodispatch_fingerprint")
  env <- new.env(parent=environment(fp))
  env$.npRmpi_training_value_projection <- function(x) list(value=x,changed=FALSE)
  incumbent <- fp; environment(incumbent) <- env
  x <- list(a=1:3,call=quote(f(x=d)))
  for (cl in list(NULL,"npregression","copula",c("user","rbandwidth"))) {
    y <- x; class(y) <- cl
    expect_identical(fp(y),incumbent(y))
  }
  expect_identical(fp(x,policy="context.v1"),incumbent(x,policy="context.v1"))
})

test_that("training terms and policy metadata remain mutation sensitive", {
  fp <- h_function(".npRmpi_autodispatch_fingerprint")
  cert <- h_function(".npRmpi_training_value_certificate")
  d <- data.frame(x=c(1,NA,3,4),y=c(2,3,4,5))
  frame <- stats::model.frame(y~x,d,na.action=stats::na.omit)
  b <- list(.np.formula.training=list(frame=frame,na.action=stats::na.omit,
                                     subset=quote(x>0)))
  before <- fp(b)
  changed <- b; attr(changed$.np.formula.training$frame,"terms") <- stats::terms(~x)
  expect_false(identical(before,fp(changed)))
  changed <- b; changed$.np.formula.training$na.action <- stats::na.exclude
  expect_false(identical(before,fp(changed)))
  changed <- b; changed$.np.formula.training$subset <- quote(x>2)
  expect_false(identical(before,fp(changed)))
  expect_false(identical(cert(list(total.time=1)),cert(list(total.time=2))))
})
