.np_entropy_integration_block_size <- function(n,
                                               target.bytes = 32 * 1024^2,
                                               max.block = 64L) {
  ## Bound the largest vectorized callback block. The conservative
  ## 64-byte-per-cell allowance covers the four distance matrices and the
  ## transient density matrices used while evaluating a block.
  max(
    1L,
    min(max.block, as.integer(floor(target.bytes / (64 * max(1, n)))))
  )
}

.np_entropy_bivariate_integral <- function(x.dat,
                                           y.dat,
                                           bw.x,
                                           bw.y,
                                           bw.joint,
                                           lower,
                                           upper,
                                           target.bytes = 32 * 1024^2,
                                           max.block = 64L) {
  n <- length(x.dat)
  block.size <- .np_entropy_integration_block_size(
    n = n,
    target.bytes = target.bytes,
    max.block = max.block
  )

  integrand <- function(xy) {
    n.eval <- ncol(xy)
    value <- numeric(n.eval)

    for (start in seq.int(1L, n.eval, by = block.size)) {
      index <- start:min(n.eval, start + block.size - 1L)
      block <- xy[, index, drop = FALSE]
      n.block <- ncol(block)

      dx <- (matrix(block[1L, ], nrow = n, ncol = n.block, byrow = TRUE) -
               x.dat) / bw.x
      dy <- (matrix(block[2L, ], nrow = n, ncol = n.block, byrow = TRUE) -
               y.dat) / bw.y
      dx.joint <- (matrix(block[1L, ], nrow = n, ncol = n.block, byrow = TRUE) -
                     x.dat) / bw.joint[1L]
      dy.joint <- (matrix(block[2L, ], nrow = n, ncol = n.block, byrow = TRUE) -
                     y.dat) / bw.joint[2L]

      f.x <- colMeans(dnorm(dx)) / bw.x
      f.y <- colMeans(dnorm(dy)) / bw.y
      f.xy <- colMeans(
        dnorm(dx.joint) * dnorm(dy.joint) /
          (bw.joint[1L] * bw.joint[2L])
      )

      value[index] <- (sqrt(f.xy) - sqrt(f.x) * sqrt(f.y))^2
    }

    matrix(value, nrow = 1L)
  }

  0.5 * adaptIntegrate(
    integrand,
    lowerLimit = lower,
    upperLimit = upper,
    vectorInterface = TRUE
  )$integral
}
