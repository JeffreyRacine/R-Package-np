test_that("fixed higher-order Gaussian convolution is an isolated fallback", {
  source_file <- test_path("..", "..", "src", "jksum.c")
  helper_file <- test_path(
    "..", "..", "src", "jksum_gaussian_fixed.c"
  )
  header_file <- test_path(
    "..", "..", "src", "jksum_gaussian_fixed.h"
  )
  skip_if_not(
    all(file.exists(c(source_file, helper_file, header_file))),
    "package C sources unavailable"
  )

  source <- paste(readLines(source_file, warn = FALSE), collapse = "\n")
  helper <- paste(readLines(helper_file, warn = FALSE), collapse = "\n")
  header <- paste(readLines(header_file, warn = FALSE), collapse = "\n")

  expect_match(
    source,
    "int fused_gaussian_fixed_higher_convolution_eligible = 0;",
    fixed = TRUE
  )
  expect_match(source, "(BANDWIDTH_reg == BW_FIXED)", fixed = TRUE)
  expect_match(source, "(operator[i] != OP_CONVOLUTION)", fixed = TRUE)
  expect_match(
    source,
    "((KERNEL_reg[i] != 1) && (KERNEL_reg[i] != 2))",
    fixed = TRUE
  )
  expect_match(
    source,
    "np_fixed_gaussian_convolution_product_try(",
    fixed = TRUE
  )
  expect_match(
    source,
    "np_convol_ckernelv(KERNEL_reg[i]",
    fixed = TRUE
  )
  expect_match(
    helper,
    "p = np_gaussian_convolution_prepare(kernel, h, hy);",
    fixed = TRUE
  )
  expect_match(
    helper,
    "np_gaussian_convolution_evaluate(&p, kernel, x-xt[i])",
    fixed = TRUE
  )
  expect_match(header, "attribute_hidden int", fixed = TRUE)

  direct <- regexpr(
    "if(fused_gaussian_convolution_kind != 0)",
    source,
    fixed = TRUE
  )
  fallback <- regexpr(
    "if(fused_gaussian_fixed_higher_convolution_eligible)",
    source,
    fixed = TRUE
  )
  expect_gt(direct, 0L)
  expect_gt(fallback, direct)
})

# Independent product-normal moment oracle, not the native Hermite expansion.
gaussian_convolution_moment_oracle <- function(x, y, hx, hy, order) {
  polynomial <- switch(as.character(order),
    "4" = c(1.5, 0, -0.5),
    "6" = c(1.875, 0, -1.25, 0, 0.125))
  variance <- hx^2 + hy^2
  mu <- (x*hy^2 + y*hx^2)/variance
  sd <- hx*hy/sqrt(variance)
  transform <- function(center, h) {
    result <- numeric(length(polynomial))
    for (k in 0:(length(polynomial)-1L))
      for (j in 0:k)
        result[j+1L] <- result[j+1L] +
          polynomial[k+1L]*choose(k,j)*((mu-center)/h)^(k-j)*(sd/h)^j
    result
  }
  a <- transform(x, hx)
  b <- transform(y, hy)
  integral <- 0
  for (i in seq_along(a)) for (j in seq_along(b)) {
    k <- i+j-2L
    moment <- if (k == 0L) 1 else if (k %% 2L) 0 else prod(seq(1,k-1,by=2))
    integral <- integral + a[i]*b[j]*moment
  }
  dnorm((x-y)/sqrt(variance))/sqrt(variance)*integral
}

test_that("Gaussian convolution and CVLS obey independent kernel calculus", {
  old <- options(np.messages=FALSE, np.largeh=FALSE, np.largelambda=FALSE,
                 np.tree=getOption("np.tree"))
  on.exit(options(old), add=TRUE)
  dat <- data.frame(a=c(-.8,-.45,-.2,.03,.22,.47,.68,.91,1.1),
                    b=c(.9,.3,-.4,.15,.7,-.6,.05,1.2,-.1))
  for (order in c(4L,6L)) for (p in 1:2) {
    x <- dat[seq_len(p)]
    h <- c(.31,.47)[seq_len(p)]
    convolution <- matrix(1,nrow(x),nrow(x))
    ordinary <- convolution
    for (d in seq_len(p)) {
      convolution <- convolution*outer(x[[d]],x[[d]],Vectorize(function(a,b)
        gaussian_convolution_moment_oracle(a,b,h[d],h[d],order)))
      z <- outer(x[[d]],x[[d]],"-")/h[d]
      polynomial <- if (order==4L) 1.5-.5*z^2 else 1.875-1.25*z^2+.125*z^4
      ordinary <- ordinary*dnorm(z)*polynomial/h[d]
    }
    diag(ordinary) <- 0
    oracle <- mean(convolution)-2*sum(ordinary)/(nrow(x)*(nrow(x)-1L))
    for (tree in c(FALSE,TRUE)) {
      options(np.tree=tree)
      actual <- npksum(txdat=x,exdat=x,bws=h,ckerorder=order,
                       operator="convolution",bandwidth.divide=TRUE)$ksum
      expect_equal(as.double(actual),rowSums(convolution),tolerance=2e-10)
      shifted <- as.data.frame(lapply(x,function(v)v+2))
      translated <- npksum(txdat=shifted,exdat=shifted,bws=h,ckerorder=order,
                          operator="convolution",bandwidth.divide=TRUE)$ksum
      expect_equal(as.double(translated),as.double(actual),tolerance=2e-10)
      scaled <- as.data.frame(lapply(x,function(v)v*3))
      rescaled <- npksum(txdat=scaled,exdat=scaled,bws=3*h,ckerorder=order,
                        operator="convolution",bandwidth.divide=TRUE)$ksum
      expect_equal(as.double(rescaled)*3^p,as.double(actual),tolerance=2e-10)
      bw <- npudensbw(dat=x,bws=h,bandwidth.compute=FALSE,
                     bwmethod="cv.ls",ckerorder=order)
      result <- npudensbw.bandwidth(dat=x,bws=bw,bandwidth.compute=TRUE,
                                   bwsolver="powell",eval.only=TRUE,nmulti=1L,
                                   powell.remin=FALSE)
      value <- if (!is.null(result$objective)) result$objective[[1L]] else -result$fval[[1L]]
      expect_equal(as.double(value),oracle,tolerance=2e-10)
    }
  }
})

test_that("adaptive Gaussian CVLS integrates the unequal-width kernels", {
  old <- options(np.messages=FALSE,np.largeh=FALSE,np.largelambda=FALSE)
  on.exit(options(old),add=TRUE)
  x <- c(-.8,-.45,-.2,.03,.22,.47,.68,.91,1.1)
  n <- length(x); k <- 3L
  h <- vapply(seq_len(n),function(i) sort(abs(x[-i]-x[i]))[k],numeric(1))
  for(order in c(4L,6L)) {
    convolution <- outer(seq_len(n),seq_len(n),Vectorize(function(i,j)
      gaussian_convolution_moment_oracle(x[i],x[j],h[i],h[j],order)))
    cross <- vapply(seq_len(n),function(i) {
      donors <- setdiff(seq_len(n),i)
      hh <- vapply(donors,function(j)
        sort(abs(x[setdiff(donors,j)]-x[j]))[k],numeric(1))
      z <- (x[i]-x[donors])/hh
      polynomial <- if(order==4L) 1.5-.5*z^2 else 1.875-1.25*z^2+.125*z^4
      mean(dnorm(z)*polynomial/hh)
    },numeric(1))
    oracle <- 2*mean(cross)-mean(convolution)
    bw <- npudensbw(dat=data.frame(x=x),bws=k,bwtype="adaptive_nn",
      ckerorder=order,bwmethod="cv.ls",bandwidth.compute=FALSE)
    actual <- npudensbw.bandwidth(dat=data.frame(x=x),bws=bw,
      bandwidth.compute=TRUE,eval.only=TRUE,nmulti=1L)$fval
    expect_equal(actual,oracle,tolerance=2e-10)
  }
})
