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
    y <- .np_hat_response(y, NULL)
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
  rows.omit.train <- .np_current_rows_omit(tdat)
  if (!is.null(y))
    rows.omit.train <- union(rows.omit.train, attr(stats::na.omit(as.data.frame(y)), "na.action"))

  if (length(rows.omit.train) > 0L)
    keep.rows[as.integer(rows.omit.train)] <- FALSE

  if (!any(keep.rows))
    stop("Training data has no rows without NAs")

  tdat <- .np_data_rows(tdat, keep.rows)
  if (!is.null(y))
    y <- y[keep.rows, , drop = FALSE]

  if (!no.e) {
    keep.eval <- rep_len(TRUE, nrow(edat))
    rows.omit.eval <- .np_current_rows_omit(edat)
    if (length(rows.omit.eval) > 0L)
      keep.eval[as.integer(rows.omit.eval)] <- FALSE
    if (!any(keep.eval))
      stop("Evaluation data has no rows without NAs")
    edat <- .np_data_rows(edat, keep.eval)
  }

  n.train <- nrow(tdat)
  target.dist <- fitted(if (no.e) {
    npudist(bws = bws, tdat = tdat)
  } else {
    npudist(bws = bws, tdat = tdat, edat = edat)
  })

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

    if (identical(bws$type, "fixed") && ncol(y) > 0L) {
      k <- ncol(y)
      sums <- .np_direct_operator_apply(
        kbw = bws, txdat = tdat, exdat = if (no.e) tdat else edat,
        operator = "integral", rhs = cbind(y, 1.0),
        where = "npudisthat direct operator apply", train.is.eval = no.e,
        separate.columns = TRUE) / n.train
      values <- sums[, seq_len(k), drop = FALSE]
      colnames(values) <- colnames(y)
      scaled <- scale_rows(values, sums[, k + 1L])
      return(if (k == 1L) as.vector(scaled) else scaled)
    }

    rhs <- if (ncol(y) == 1L) {
      y
    } else {
      lapply(seq_len(ncol(y)), function(j) y[, j, drop = FALSE])
    }

    fixed <- identical(bws$type, "fixed")
    operator.apply <- if (fixed) .np_direct_operator_apply else .np_exact_operator_apply
    label <- if (fixed) "npudisthat direct operator" else "npudisthat exact operator"
    apply.column <- function(column, operation) {
      operator.apply(kbw = bws, txdat = tdat, exdat = if (no.e) tdat else edat,
                     operator = "integral", rhs = column,
                     where = paste(label, operation), train.is.eval = no.e) / n.train
    }

    # All response columns have the same kernel row normalization.
    probe.sum <- apply.column(matrix(1.0, n.train, 1L), "probe")
    apply.one <- function(column) {
      out <- apply.column(column, "apply")
      scale_rows(out, probe.sum)
    }
    out <- if (is.list(rhs)) do.call(cbind, lapply(rhs, apply.one)) else apply.one(rhs)

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
      where = "npudisthat direct operator",
      train.is.eval = no.e
    ) / n.train
  } else {
    H <- .np_exact_operator_matrix(
      kbw = bws,
      txdat = tdat,
      exdat = if (no.e) tdat else edat,
      operator = "integral",
      where = "npudisthat exact operator",
      train.is.eval = no.e
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
