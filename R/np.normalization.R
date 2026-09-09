# A normalizing sum is not an additive estimate: never replace it by epsilon.
# zero.rows is lazy evidence used only for an exceptional external point row.
.np_normalization_denominator <- function(denominator, where,
                                          allow.empty.rows = FALSE,
                                          zero.rows = NULL,
                                          replicate.offset = NULL) {
  bad <- which(!is.finite(denominator) | denominator == 0.0)
  if (!length(bad))
    return(denominator)

  empty <- integer()
  if (isTRUE(allow.empty.rows) && is.null(dim(denominator))) {
    zero <- bad[is.finite(denominator[bad]) & denominator[bad] == 0.0]
    if (length(zero)) {
      evidence <- zero.rows
      if (!is.null(evidence))
        empty <- zero[which(evidence[zero] %in% TRUE)]
    }
  }
  invalid <- setdiff(bad, empty)
  if (length(invalid)) {
    first <- invalid[[1L]]
    if (is.matrix(denominator)) {
      at <- arrayInd(first, dim(denominator))
      position <- if (is.null(replicate.offset)) {
        sprintf("row %d, column %d", at[1L], at[2L])
      } else {
        sprintf("bootstrap replication %d, evaluation row %d",
                replicate.offset + at[1L], at[2L])
      }
    } else {
      position <- sprintf("evaluation row %d", first)
      if (!is.null(replicate.offset))
        position <- sprintf("bootstrap replication %d, %s",
                            replicate.offset + 1L, position)
    }
    reason <- if (!is.finite(denominator[first])) "non-finite" else "zero"
    stop(sprintf("%s: %s normalizing weight sum at %s; the required ratio is undefined",
                 where, reason, position), call. = FALSE)
  }
  denominator[empty] <- NA_real_
  denominator
}

.np_bootstrap_ratio <- function(numerator, denominator, where,
                                 first.replication = 1L) {
  numerator / .np_normalization_denominator(
    denominator, where, replicate.offset =
      if (is.null(first.replication)) NULL else first.replication - 1L)
}

.np_normalization_finish <- function(value, denominator, where,
                                     defer.empty.rows = FALSE) {
  if (!anyNA(denominator))
    return(value)
  empty <- as.integer(is.na(denominator))
  if (isTRUE(defer.empty.rows)) {
    attr(value, ".np.empty.rows") <- empty
  } else {
    rows <- which(empty == 1L)
    .np_warning(sprintf(
      "%s: all computed kernel weights are zero at %d external evaluation row(s) (%s%s); returning NA for those rows",
      where, length(rows), paste(utils::head(rows, 8L), collapse = ", "),
      if (length(rows) > 8L) ", ..." else ""), call. = FALSE)
  }
  value
}
