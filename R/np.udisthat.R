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
  target.dist <- fitted(npudist(
    bws = bws,
    tdat = tdat,
    edat = if (no.e) tdat else edat
  ))

  scale_rows <- function(x, probe.sum) {
    x <- as.matrix(x)
    probe.sum <- as.vector(probe.sum)
    safe.sum <- probe.sum
    safe.sum[abs(safe.sum) < .Machine$double.xmin] <- .Machine$double.xmin
    scaled <- sweep(x, 1L, target.dist / safe.sum, "*")
    if (ncol(scaled) == 1L)
      scaled[, 1L, drop = FALSE]
    else
      scaled
  }

  if (identical(output, "apply")) {
    if (is.null(y))
      stop("argument 'y' is required when output='apply'")

    rhs <- if (ncol(y) == 1L) {
      y
    } else {
      lapply(seq_len(ncol(y)), function(j) y[, j, drop = FALSE])
    }

    apply_one <- function(rhs.col) {
      if (identical(bws$type, "fixed")) {
        out <- .np_direct_operator_apply(
          kbw = bws,
          txdat = tdat,
          exdat = if (no.e) tdat else edat,
          operator = "integral",
          rhs = rhs.col,
          where = "npudisthat direct operator apply"
        ) / n.train
        probe.sum <- .np_direct_operator_apply(
          kbw = bws,
          txdat = tdat,
          exdat = if (no.e) tdat else edat,
          operator = "integral",
          rhs = matrix(1.0, nrow = n.train, ncol = 1L),
          where = "npudisthat direct operator probe"
        ) / n.train
      } else {
        out <- .np_exact_operator_apply(
          kbw = bws,
          txdat = tdat,
          exdat = if (no.e) tdat else edat,
          operator = "integral",
          rhs = rhs.col,
          where = "npudisthat exact operator apply"
        ) / n.train
        probe.sum <- .np_exact_operator_apply(
          kbw = bws,
          txdat = tdat,
          exdat = if (no.e) tdat else edat,
          operator = "integral",
          rhs = matrix(1.0, nrow = n.train, ncol = 1L),
          where = "npudisthat exact operator probe"
        ) / n.train
      }

      scale_rows(out, probe.sum)
    }

    out <- if (is.list(rhs)) {
      do.call(cbind, lapply(rhs, apply_one))
    } else {
      apply_one(rhs)
    }

    if (ncol(out) == 1L)
      return(as.vector(out))
    return(out)
  }

  if (identical(bws$type, "fixed")) {
    H <- .np_direct_operator_matrix(
      kbw = bws,
      txdat = tdat,
      exdat = if (no.e) tdat else edat,
      operator = "integral",
      where = "npudisthat direct operator"
    ) / n.train
  } else {
    H <- .np_exact_operator_matrix(
      kbw = bws,
      txdat = tdat,
      exdat = if (no.e) tdat else edat,
      operator = "integral",
      where = "npudisthat exact operator"
    ) / n.train
  }

  H <- scale_rows(H, rowSums(H))

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
