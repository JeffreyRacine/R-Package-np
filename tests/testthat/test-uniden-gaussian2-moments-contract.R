test_that("Gaussian2 coefficients retain the continuous moment-ratio limit", {
  owner <- npuniden.boundary
  parts <- as.list(body(owner))
  location <- which(vapply(parts,function(z)is.call(z) &&
    identical(z[[1]],as.name("<-")) && identical(z[[2]],as.name("int.kernel.squared")),FALSE))
  body(owner) <- as.call(c(parts[seq_len(location-1L)],list(quote(return(environment())))))
  X <- seq(.03,.96,length.out=20)^1.3
  kernel <- owner(X,Y=.5,h=.18,a=0,b=1,kertype="gaussian2")$kernel
  # Independent integration of the original truncated-normal moments.
  for(h in c(.18,.5,2,10,243,1000,9000)) for(y in c(.07,.31,.67,.91)) {
    lower <- -y/h; upper <- (1-y)/h
    m <- vapply(0:3,function(k)integrate(function(z)z^k*dnorm(z),lower,upper,
      rel.tol=1e-12,abs.tol=1e-13)$value,0.0)
    mu <- m/m[1L];den <- mu[4L]-mu[2L]*mu[3L]
    z <- (y-X)/h
    oracle <- (mu[4L]-mu[2L]*z^2)/den*dnorm(z)/(h*m[1L])
    expect_equal(kernel(y,X,h,0,1),oracle,tolerance=1e-10)
  }
  for(h in c(.005,.04,.18,.5,2)) {
    c <- .5/h
    mass <- pnorm(c)-pnorm(-c)
    m2 <- integrate(function(z)z^2*dnorm(z),-c,c,rel.tol=1e-12)$value/mass
    z <- (.5-X)/h
    oracle <- (c^2-z^2)/(c^2-m2)*dnorm(z)/(h*mass)
    expect_equal(kernel(.5,X,h,0,1),oracle,tolerance=1e-12)
    for(direction in c(-1,1)) {
      y <- .5+direction*h*1e-12/max(1,c)
      expect_equal(kernel(y,X,h,0,1),oracle,tolerance=1e-9)
    }
    for(width in c(.05,4))
      expect_equal(width*kernel(7+width*.5,7+width*X,h*width,7,7+width),
                   oracle,tolerance=1e-9)
  }
  # Both sides of the series/expm1 seam must agree continuously.
  h <- .5
  for(sign in c(-1,1)) {
    y <- .5+sign*h*1e-3/(.5/h)
    below <- kernel(.5+(y-.5)*(1-1e-8),X,h,0,1)
    above <- kernel(.5+(y-.5)*(1+1e-8),X,h,0,1)
    expect_equal(below,above,tolerance=1e-9)
  }
  for(proper in c(FALSE,TRUE)) {
    fit <- npuniden.boundary(X,Y=c(.25,.5,.75),h=.18,a=0,b=1,
                            kertype="gaussian2",proper=proper)
    expect_true(all(is.finite(unlist(fit[c("f","F","sd.f","sd.F")]))))
  }
  for(method in c("cv.ls","cv.ml")) {
    fit <- npuniden.boundary(c(X,.5),Y=c(.25,.5,.75),a=0,b=1,
      grid=c(.08,.18,.4),bwmethod=method,kertype="gaussian2")
    expect_true(all(is.finite(unlist(fit))))
  }
})

test_that("Gaussian2 preserves its moment kernel through the large-bandwidth limit", {
  owner <- npuniden.boundary
  parts <- as.list(body(owner))
  location <- which(vapply(parts, function(z) is.call(z) &&
    identical(z[[1L]], as.name("<-")) &&
    identical(z[[2L]], as.name("int.kernel.squared")), FALSE))
  body(owner) <- as.call(c(parts[seq_len(location - 1L)], list(quote(return(environment())))))
  X <- seq(.03, .96, length.out = 20)^1.3
  kernel <- owner(X, Y = .5, h = 1, a = 0, b = 1, kertype = "gaussian2")$kernel
  for (h in c(9.999, 10.001, 9999, 10000, 10001, 1e8, 1e100, 1e300))
    for (x in c(0, .07, .5, .91, 1)) {
      # Independent numerical moments on physical unit support avoid the same
      # tiny-standardized-moment algebra as the implementation.
      moment <- vapply(0:3, function(k) integrate(function(w)
        w^k * exp(-.5 * (w / h)^2), -x, 1 - x,
        rel.tol = 1e-12, abs.tol = 1e-13)$value, numeric(1L))
      ratio <- if (x == .5) .25 else moment[4L] / moment[2L]
      oracle <- (ratio - (x - X)^2) / (ratio * moment[1L] - moment[3L]) *
        exp(-.5 * ((x - X) / h)^2)
      expect_equal(kernel(x, X, h, 0, 1), oracle, tolerance = 1e-11)
      expect_equal(integrate(function(z) kernel(x, z, h, 0, 1), 0, 1,
        rel.tol = 1e-11)$value, 1, tolerance = 1e-11)
      expect_equal(integrate(function(z) (z - x) * kernel(x, z, h, 0, 1), 0, 1,
        abs.tol = 1e-12)$value, 0, tolerance = 1e-11)
      expect_equal(4 * kernel(7 + 4*x, 7 + 4*X, 4*h, 7, 11), oracle,
        tolerance = 1e-11)
    }
  for (x in c(0, .07, .5, .91, 1)) {
    limit <- 6 * (.25 + (x - .5)^2 - (x - X)^2)
    expect_equal(kernel(x, X, 1e300, 0, 1), limit, tolerance = 1e-13)
    expect_equal(kernel(x, X, 10000*(1 - 1e-8), 0, 1),
      kernel(x, X, 10000*(1 + 1e-8), 0, 1), tolerance = 1e-11)
  }
})

test_that("Gaussian2 CV objectives retain the moment estimator beyond the old seam", {
  owner <- npuniden.boundary
  parts <- as.list(body(owner))
  location <- which(vapply(parts, function(z) is.call(z) &&
    identical(z[[1L]], as.name("if")) &&
    identical(z[[2L]], quote(is.null(h) && cv == "grid-hybrid")), FALSE))
  body(owner) <- as.call(c(parts[seq_len(location - 1L)], list(quote(return(environment())))))
  X <- seq(.03, .96, length.out = 20)^1.3
  ns <- environment(npuniden.boundary)
  total <- get(".np_quadrature_total", ns)
  oracle.kernel <- function(x, h) {
    m <- vapply(0:3, function(k) integrate(function(w) w^k * exp(-.5*(w/h)^2),
      -x, 1-x, rel.tol = 1e-12, abs.tol = 1e-13)$value, numeric(1L))
    ratio <- m[4L]/m[2L]
    (ratio-(x-X)^2)/(ratio*m[1L]-m[3L]) * exp(-.5*((x-X)/h)^2)
  }
  for (method in c("cv.ls", "cv.ml")) {
    context <- owner(X, a = 0, b = 1, kertype = "gaussian2", bwmethod = method)
    values <- numeric(3L); index <- 0L
    for (h in c(9999, 10001, 1e8)) {
      weights <- vapply(X, oracle.kernel, numeric(length(X)), h = h)
      loo <- (colSums(weights)-diag(weights))/(length(X)-1L)
      oracle <- if (method == "cv.ml") sum(log(pmax(loo, .Machine$double.xmin))) else
        total(X, colMeans(weights)^2)-2*mean(loo)
      actual <- context$cv.function(h, X, 0, 1)
      expect_equal(actual, oracle, tolerance = 1e-11)
      index <- index+1L; values[index] <- actual
    }
    expect_lt(abs(values[1L]-values[2L]), 1e-8)
  }
})
