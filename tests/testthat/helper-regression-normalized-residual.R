# Independent ordinary-regression variance oracle for the retained training map.
# Do not use the native residual helper here: explicit I - H also covers an
# actual ridge map that does not reproduce constants exactly.
hc0_normalized_training_residual <- function(training.hat, response,
                                             training.mean = NULL,
                                             constant.reproduction = FALSE) {
  response <- as.double(response)
  training.hat <- as.matrix(training.hat)
  n <- length(response)
  stopifnot(identical(dim(training.hat), c(n, n)),
            all(is.finite(training.hat)), all(is.finite(response)))
  if (isTRUE(constant.reproduction)) {
    # The caller must establish this owner contract independently. Complete
    # the residual diagonal from the actual off-diagonal map, not 1-H_ii,
    # and center before multiplying so near-interpolation remains informative.
    return(vapply(seq_len(n), function(i) {
      weights <- training.hat[i, -i]
      scale <- max(abs(weights))
      if (scale == 0)
        return(NA_real_)
      weights <- weights / scale
      sum(weights * (response[i] - response[-i])) /
        sqrt(sum(weights)^2 + sum(weights^2))
    }, numeric(1L)))
  }
  fitted <- if (is.null(training.mean)) {
    drop(training.hat %*% response)
  } else {
    stopifnot(length(training.mean) == n, all(is.finite(training.mean)))
    as.double(training.mean)
  }
  residual.map <- diag(n) - training.hat
  row.norm <- apply(residual.map, 1L, function(row) {
    scale <- max(abs(row))
    if (scale == 0)
      return(0)
    scale * sqrt(sum((row / scale)^2))
  })
  normalized <- rep.int(NA_real_, n)
  identified <- row.norm > 0
  normalized[identified] <- (response - fitted)[identified] /
    row.norm[identified]
  # q == 0 means unavailable donor variance, never an epsilon replacement.
  # Direction-specific zero/dependency contracts are tested separately.
  normalized
}
