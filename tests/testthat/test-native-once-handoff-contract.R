native_once_fixture <- function() {
  set.seed(7191)
  x <- data.frame(x = runif(32), u = factor(rep(0:1, 16)))
  list(x = x, y = x$x + as.numeric(x$u) + rnorm(32, sd = .3),
       z = data.frame(z = runif(32)),
       xd = data.frame(x = x$x, o = ordered(rep(0:1, 16))))
}

test_that("native one-call selectors consume the original data promises once", {
  f <- native_once_fixture()
  x <- f$x; y <- f$y; z <- f$z; xd <- f$xd
  xc <- x[1L]; dfy <- data.frame(y = y)
  cases <- list(
    reg = quote(npreg(c(.4, .2), touch("x", x), touch("y", y))),
    density = quote(npudens(c(.4, .2), touch("x", x))),
    distribution = quote(npudist(bws = c(.4, .2), tdat = touch("x", xd))),
    cdensity = quote(npcdens(bws = c(.4, .2, .4), txdat = touch("x", x),
                             tydat = touch("y", dfy))),
    cdistribution = quote(npcdist(bws = c(.4, .2, .4), txdat = touch("x", x),
                                   tydat = touch("y", dfy))),
    index = quote(npindex(bws = c(.4, .5, 1), txdat = touch("x", x),
                          tydat = touch("y", y), se = FALSE)),
    plreg = quote(npplreg(bws = matrix(.4, 3, 1), txdat = touch("x", x),
                         tydat = touch("y", y), tzdat = touch("z", z))),
    scoef = quote(npscoef(bws = .4, txdat = touch("x", xc),
                         tydat = touch("y", y), tzdat = touch("z", z))),
    qreg = quote(npqreg(bws = c(.4, .2, .4), txdat = touch("x", x),
                       tydat = touch("y", dfy))),
    conmode = quote(npconmode(bws = c(.4, .2, .2), txdat = touch("x", x),
                             tydat = touch("y", data.frame(y = x$u)))),
    sigtest = quote(npsigtest(bws = c(.4, .2), xdat = touch("x", x),
                             ydat = touch("y", y), B = 9, random.seed = 3)),
    lsqreg = quote(nplsqreg(bws = c(.4, .2), txdat = touch("x", x),
                           tydat = touch("y", y), bandwidth.compute = FALSE)))
  strip <- function(e) {
    if (!is.call(e)) return(e)
    if (identical(e[[1]], as.name("touch"))) return(strip(e[[3]]))
    for (i in seq_along(e)[-1L]) e[[i]] <- strip(e[[i]])
    e
  }
  fields <- function(a) {
    b <- a[["bws", exact = TRUE]]
    list(point = if (!is.null(a[["P", exact = TRUE]])) a$P else fitted(a),
         statistic = a[["In", exact = TRUE]],
         names = b[intersect(c("xnames", "ynames", "znames", "names"), names(b))],
         training = b[[".np.native.training", exact = TRUE]])
  }
  for (id in names(cases)) for (wrapped in c(FALSE, TRUE)) {
    counts <- c(x = 0L, y = 0L, z = 0L)
    touch <- function(role, value) {
      counts[role] <<- counts[role] + 1L
      invisible(runif(1))
      if (counts[role] > 1L) {
        if (is.data.frame(value)) value[[1L]] <- rev(value[[1L]]) else
          value <- rev(value)
      }
      value
    }
    expr <- cases[[id]]
    if (wrapped) {
      w <- function(...) NULL
      body(w) <- as.call(list(expr[[1L]], as.name("...")))
      expr[[1L]] <- as.name("w")
    }
    set.seed(2908)
    actual <- eval(expr)
    actual.rng <- .Random.seed
    expect_true(all(counts %in% 0:1), info = paste(id, wrapped))
    set.seed(2908)
    invisible(runif(sum(counts > 0L)))
    expected <- eval(strip(expr))
    expect_identical(fields(actual), fields(expected), info = paste(id, wrapped))
    expect_identical(.Random.seed, actual.rng, info = paste(id, wrapped, "RNG"))
  }
})

test_that("native errors and warnings do not replay data expressions", {
  f <- native_once_fixture()
  count <- 0L
  bad <- function() { count <<- count + 1L; stop("native-once-sentinel") }
  expect_error(npreg(bws = c(.4, .2), txdat = bad(), tydat = f$y),
               "native-once-sentinel")
  expect_identical(count, 1L)
  count <- 0L
  warn <- function() { count <<- count + 1L; warning("native-once-warning"); f$x }
  expect_warning(npreg(bws = c(.4, .2), txdat = warn(), tydat = f$y),
                 "native-once-warning")
  expect_identical(count, 1L)
})
