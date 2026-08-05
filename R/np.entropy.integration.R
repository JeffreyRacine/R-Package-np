.np_entropy_bivariate_integral <- function(x.dat,
                                           y.dat,
                                           bw.x,
                                           bw.y,
                                           bw.joint,
                                           lower,
                                           upper) {
  x.dat <- as.double(x.dat)
  y.dat <- as.double(y.dat)
  bandwidths <- as.double(c(bw.x, bw.y, bw.joint))

  integrand <- function(xy) {
    .Call(
      "C_np_entropy_gaussian_integrand",
      xy,
      x.dat,
      y.dat,
      bandwidths,
      PACKAGE = "np"
    )
  }

  0.5 * adaptIntegrate(
    integrand,
    lowerLimit = lower,
    upperLimit = upper,
    vectorInterface = TRUE
  )$integral
}
