.npcdenshat_make_xbw <- function(bws, txdat) {
  xbw.args <- list(
    xdat = txdat,
    ydat = rep.int(0.0, nrow(txdat)),
    bws = bws$xbw,
    regtype = bws$regtype,
    bwtype = bws$type,
    bandwidth.compute = FALSE,
    ckertype = bws$cxkertype,
    ckerorder = bws$cxkerorder,
    ckerbound = bws$cxkerbound,
    ckerlb = bws$cxkerlb,
    ckerub = bws$cxkerub,
    ukertype = bws$uxkertype,
    okertype = bws$oxkertype
  )

  if (!is.null(bws$basis))
    xbw.args$basis <- bws$basis
  if (!is.null(bws$degree))
    xbw.args$degree <- bws$degree
  if (!is.null(bws$bernstein.basis))
    xbw.args$bernstein.basis <- bws$bernstein.basis

  do.call(npregbw, xbw.args)
}

.npcdenshat_make_ybw <- function(bws, tydat) {
  kbandwidth.numeric(
    bw = bws$ybw,
    bwscaling = FALSE,
    bwtype = bws$type,
    ckertype = bws$cykertype,
    ckerorder = bws$cykerorder,
    ckerbound = bws$cykerbound,
    ckerlb = bws$cykerlb,
    ckerub = bws$cykerub,
    ukertype = bws$uykertype,
    okertype = bws$oykertype,
    nobs = nrow(tydat),
    xdati = untangle(tydat),
    ydati = NULL,
    xnames = names(tydat),
    ynames = NULL
  )
}

.npcdenshat_make_xhat <- function(xbw, txdat, exdat) {
  n.train <- nrow(txdat)

  Hx <- vapply(
    seq_len(n.train),
    function(j) {
      npreg(
        bws = xbw,
        txdat = txdat,
        tydat = as.numeric(seq_len(n.train) == j),
        exdat = exdat,
        warn.glp.gradient = FALSE
      )$mean
    },
    numeric(nrow(exdat))
  )

  Hx <- as.matrix(Hx)
  if (nrow(Hx) != nrow(exdat))
    Hx <- t(Hx)
  Hx
}

npcdenshat <- function(bws,
                       txdat = stop("training data 'txdat' missing"),
                       tydat = stop("training data 'tydat' missing"),
                       exdat,
                       eydat,
                       y = NULL,
                       output = c("matrix", "apply")) {
  output <- match.arg(output)

  if (!inherits(bws, "conbandwidth")) {
    stop("argument 'bws' must inherit from class 'conbandwidth' in npcdenshat()")
  }
  if (!identical(bws$type, "fixed")) {
    stop("npcdenshat() currently supports bwtype = 'fixed' only")
  }

  if (xor(missing(exdat), missing(eydat)))
    stop("evaluation data must be supplied for both 'exdat' and 'eydat'")

  no.exy <- missing(exdat)

  txdat <- toFrame(txdat)
  tydat <- toFrame(tydat)

  if (!no.exy) {
    exdat <- toFrame(exdat)
    eydat <- toFrame(eydat)

    if (!(txdat %~% exdat))
      stop("'txdat' and 'exdat' are not similar data frames!")
    if (!(tydat %~% eydat))
      stop("'tydat' and 'eydat' are not similar data frames!")
  }

  if (!is.null(y)) {
    if (is.factor(y) || is.vector(y)) {
      y <- matrix(as.double(y), ncol = 1L)
    } else {
      y <- as.matrix(y)
      storage.mode(y) <- "double"
    }

    if (nrow(y) != nrow(txdat))
      stop("number of rows in 'y' must match the number of training rows in 'txdat'")
  }

  keep.rows <- rep_len(TRUE, nrow(txdat))
  rows.omit.train <- attr(stats::na.omit(data.frame(txdat, tydat)), "na.action")
  if (!is.null(y))
    rows.omit.train <- union(rows.omit.train, attr(stats::na.omit(as.data.frame(y)), "na.action"))

  if (length(rows.omit.train) > 0L)
    keep.rows[as.integer(rows.omit.train)] <- FALSE

  if (!any(keep.rows))
    stop("Training data has no rows without NAs")

  txdat <- txdat[keep.rows, , drop = FALSE]
  tydat <- tydat[keep.rows, , drop = FALSE]
  if (!is.null(y))
    y <- y[keep.rows, , drop = FALSE]

  if (!no.exy) {
    keep.eval <- rep_len(TRUE, nrow(exdat))
    rows.omit.eval <- attr(stats::na.omit(data.frame(exdat, eydat)), "na.action")
    if (length(rows.omit.eval) > 0L)
      keep.eval[as.integer(rows.omit.eval)] <- FALSE

    if (!any(keep.eval))
      stop("Evaluation data has no rows without NAs")

    exdat <- exdat[keep.eval, , drop = FALSE]
    eydat <- eydat[keep.eval, , drop = FALSE]
  } else {
    exdat <- txdat
    eydat <- tydat
  }

  xbw <- .npcdenshat_make_xbw(bws = bws, txdat = txdat)
  ybw <- .npcdenshat_make_ybw(bws = bws, tydat = tydat)

  Hx <- .npcdenshat_make_xhat(xbw = xbw, txdat = txdat, exdat = exdat)
  Gy <- npksum(
    bws = ybw,
    txdat = tydat,
    exdat = eydat,
    tydat = diag(nrow(txdat)),
    operator = "normal",
    bandwidth.divide = TRUE
  )$ksum

  H <- as.matrix(Hx) * as.matrix(Gy)

  if (identical(output, "apply")) {
    if (is.null(y))
      stop("argument 'y' is required when output='apply'")

    out <- H %*% y
    if (ncol(out) == 1L)
      return(as.vector(out))
    return(out)
  }

  class(H) <- c("npcdenshat", "matrix")
  attr(H, "bws") <- bws
  attr(H, "txdat") <- txdat
  attr(H, "tydat") <- tydat
  attr(H, "exdat") <- exdat
  attr(H, "eydat") <- eydat
  attr(H, "trainiseval") <- no.exy
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
