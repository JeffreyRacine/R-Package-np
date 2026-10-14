# Tiny native numerical contracts; public/MPI/kernel matrices stay in focused
# validation. Existing paired-SE tests cover demand and ordinary route behavior.
pair_center_native <- function(z, d, w, r, slope=NULL) {
  .Call("C_np_conditional_lp_pair_se", z, d, w, r, slope, PACKAGE="np")
}

pair_center_mass_oracle <- function(z, d, w, r, ridge=FALSE) {
  n <- nrow(z)
  endpoint <- function(mass, j) {
    a <- crossprod(z, z*(mass*w[,j]))
    b <- crossprod(z, mass*w[,j]*r[,j])
    rho <- if(ridge)max(abs(diag(a)))/n else 0
    if(ridge)b[1L] <- b[1L]*(1+rho/a[1L,1L])
    drop(d[j,] %*% solve(a+diag(rho,ncol(z)),b))
  }
  eps <- 1e-5
  influence <- vapply(seq_len(n),function(i) {
    up <- down <- rep(1,n);up[i] <- 1+eps;down[i] <- 1-eps
    (endpoint(up,1L)-endpoint(up,2L)-endpoint(down,1L)+endpoint(down,2L))/(2*eps)
  },0.)
  list(point=vapply(1:2,function(j)endpoint(rep(1,n),j),0.),
    se=sqrt(n/(n-1)*sum((influence-mean(influence))^2)))
}

test_that("paired conditional constants use exact endpoint response translations", {
  x <- seq(-.4,.4,length.out=8L)
  w <- cbind(c(rep(1.,4L),rep(0.,4L)),c(rep(0.,4L),rep(1.,4L)))
  for(ridged in c(FALSE,TRUE)) {
    z <- if(ridged)cbind(1,x,x) else cbind(1,x)
    d <- matrix(c(1,1,rep(0,2*(ncol(z)-1L))),2L)
    for(value in c(0,1,.3)) {
      r <- cbind(rep(value,8L),rep(.8,8L))
      # Unsupported donors do not define an endpoint's response anchor.
      r[w==0] <- -7
      got <- pair_center_native(z,d,w,r)
      expect_identical(got[1L,6L],0)
      expect_identical(got[1L,1L],0)
      expect_identical(got[1L,2:3],c(value,.8))
    }
  }
  # Restoration uses the actual evaluation intercept, including derivative rows.
  z <- cbind(1,x);d <- rbind(c(2,0),c(0,1))
  got <- pair_center_native(z,d,w,matrix(.3,8L,2L))
  expect_identical(got[1L,c(1L,2L,3L,6L)],c(0,.6,0,0))
})

test_that("centered paired scores retain all accepted ridge mass derivatives", {
  x <- seq(-.4,.4,length.out=8L)
  w <- cbind(seq(.4,1.1,length.out=8L),rev(seq(.3,1.,length.out=8L)))
  r <- cbind(.3+sin(3*x),.7+cos(2*x))
  for(ridged in c(FALSE,TRUE)) {
    z <- if(ridged)cbind(1,x,x) else cbind(1,x)
    d <- if(ridged)rbind(c(1,.1,.1),c(1,-.15,-.15)) else rbind(c(1,.1),c(1,-.15))
    oracle <- pair_center_mass_oracle(z,d,w,r,ridged)
    got <- pair_center_native(z,d,w,r)
    expect_identical(got[1L,6L],0)
    expect_equal(got[1L,2:3],oracle$point,tolerance=1e-12)
    expect_equal(got[1L,1L],oracle$se,tolerance=1e-8)
    shifted <- pair_center_native(z,d,w,sweep(r,2,c(16,-8),'+'))
    expect_equal(shifted[1L,2:3],got[1L,2:3]+c(16,-8),tolerance=1e-12)
    expect_equal(shifted[1L,1L],got[1L,1L],tolerance=1e-12)
    scaled <- pair_center_native(z,d,w,-2*r)
    expect_equal(scaled[1L,2:3],-2*got[1L,2:3],tolerance=1e-12)
    expect_equal(scaled[1L,1L],2*got[1L,1L],tolerance=1e-12)
    ix <- c(8L,1L,7L,2L,6L,3L,5L,4L)
    permuted <- pair_center_native(z[ix,,drop=FALSE],d,w[ix,,drop=FALSE],r[ix,,drop=FALSE])
    expect_equal(permuted,got,tolerance=1e-12)
  }
})

test_that("paired centering preserves nonunit and finite-range exclusions", {
  z <- matrix(c(1.,2.,3.),3L);d <- matrix(1.,2L)
  got <- pair_center_native(z,d,matrix(1.,3L,2L),cbind(rep(2.,3L),rep(0.,3L)))
  expect_identical(got[1L,6L],0)
  expect_equal(got[1L,2L],6/7,tolerance=1e-14)
  expect_equal(got[1L,1L],sqrt(84)/49,tolerance=1e-14)
  # The raw finite map remains available when subtracting its endpoint anchor
  # would overflow. No tolerant-zero or large-value clamp is permitted.
  # Powers of two make raw weighted products exactly +/-1 on every IEEE
  # BLAS, while the attempted response difference still overflows.
  huge <- 2^1023
  got <- pair_center_native(matrix(1.,2L),matrix(1.,2L),matrix(2^-1023,2L,2L),
    cbind(c(huge,-huge),c(0.,0.)))
  expect_identical(got[1L,6L],0)
  expect_equal(got[1L,1L]/huge,1,tolerance=2e-12)
  expect_identical(got[1L,2:3],c(0,0))
  r <- cbind(.3+1e-8*c(1,-2,3,-4),rep(.7,4L))
  got <- pair_center_native(matrix(1.,4L),matrix(1.,2L),matrix(1.,4L,2L),r)
  expect_identical(got[1L,6L],0)
  expect_gt(got[1L,1L],0)
  expect_equal(got[1L,1L],sd(r[,1L])/2,tolerance=1e-12)
})

test_that("constant paired endpoints retain refusal and recovery contracts", {
  z <- cbind(1,seq(-.4,.4,length.out=8L));d <- rbind(c(1,0),c(1,.1))
  w <- matrix(1.,8L,2L);r <- matrix(.3,8L,2L)
  valid <- pair_center_native(z,d,w,r)
  empty <- pair_center_native(z,d,replace(w,seq_len(8L),0),r)
  expect_identical(empty[1L,6L],1)
  expect_true(is.na(empty[1L,1L]))
  expect_error(pair_center_native(z,d,w,replace(r,1L,NA_real_)),"non-finite input")
  expect_identical(pair_center_native(z,d,w,r),valid)
  density <- pair_center_native(z,d,w,r,matrix(0.,8L,2L))
  expect_identical(density[1L,6L],4)
  expect_true(is.na(density[1L,1L]))
  # Forced ridge, exact maximum-scale tie with unequal mass derivatives.
  zn <- cbind(1,c(0.,2.),c(0.,2.));dn <- rbind(c(1,0,0),c(1,0,0))
  tied <- pair_center_native(zn,dn,matrix(c(3.,1.),2L,2L),matrix(.3,2L,2L))
  expect_identical(tied[1L,6L],2)
  expect_true(is.na(tied[1L,1L]))
  expect_identical(pair_center_native(z,d,w,r),valid)
})
