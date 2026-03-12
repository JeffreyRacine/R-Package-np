npudisthat <- function(bws,
                       tdat = stop("training data 'tdat' missing"),
                       edat,
                       y = NULL,
                       output = c("matrix", "apply")) {
  output <- match.arg(output)

  if (!inherits(bws, "dbandwidth")) {
    stop("argument 'bws' must inherit from class 'dbandwidth' in npudisthat()")
  }

  tdat <- toFrame(tdat)
  no.e <- missing(edat)

  if (!no.e) {
    edat <- toFrame(edat)
    if (!(tdat %~% edat))
      stop("'tdat' and 'edat' are not similar data frames!")
  }

  if (!is.null(y)) {
    if (is.factor(y) || is.vector(y)) {
      y <- matrix(as.double(y), ncol = 1L)
    } else {
      y <- as.matrix(y)
      storage.mode(y) <- "double"
    }

    if (nrow(y) != nrow(tdat))
      stop("number of rows in 'y' must match the number of training rows in 'tdat'")
  }

  keep.rows <- rep_len(TRUE, nrow(tdat))
  rows.omit.train <- attr(stats::na.omit(tdat), "na.action")
  if (!is.null(y))
    rows.omit.train <- union(rows.omit.train, attr(stats::na.omit(as.data.frame(y)), "na.action"))

  if (length(rows.omit.train) > 0L)
    keep.rows[as.integer(rows.omit.train)] <- FALSE

  if (!any(keep.rows))
    stop("Training data has no rows without NAs")

  tdat <- tdat[keep.rows, , drop = FALSE]
  if (!is.null(y))
    y <- y[keep.rows, , drop = FALSE]

  if (!no.e) {
    keep.eval <- rep_len(TRUE, nrow(edat))
    rows.omit.eval <- attr(stats::na.omit(edat), "na.action")
    if (length(rows.omit.eval) > 0L)
      keep.eval[as.integer(rows.omit.eval)] <- FALSE
    if (!any(keep.eval))
      stop("Evaluation data has no rows without NAs")
    edat <- edat[keep.eval, , drop = FALSE]
  }

  n.train <- nrow(tdat)
  if (identical(bws$type, "fixed")) {
    H <- .np_direct_operator_matrix(
      kbw = bws,
      txdat = tdat,
      exdat = if (no.e) tdat else edat,
      operator = "integral",
      where = "npudisthat direct operator"
    ) / n.train
  } else {
    H <- npksum(
      bws = bws,
      txdat = tdat,
      exdat = if (no.e) tdat else edat,
      tydat = diag(n.train),
      operator = "integral",
      bandwidth.divide = TRUE
    )$ksum / n.train
    H <- as.matrix(H)
  }

  if (identical(output, "apply")) {
    if (is.null(y))
      stop("argument 'y' is required when output='apply'")

    out <- H %*% y
    if (ncol(out) == 1L)
      return(as.vector(out))
    return(out)
  }

  class(H) <- c("npudisthat", "matrix")
  attr(H, "bws") <- bws
  attr(H, "tdat") <- tdat
  attr(H, "edat") <- if (no.e) tdat else edat
  attr(H, "trainiseval") <- no.e
  attr(H, "rows.omit.train") <- rows.omit.train
  attr(H, "call") <- match.call(expand.dots = FALSE)

  if (!is.null(y)) {
    Hy <- H %*% y
    if (ncol(Hy) == 1L)
      Hy <- as.vector(Hy)
    attr(H, "Hy") <- Hy
  }

  H
}
