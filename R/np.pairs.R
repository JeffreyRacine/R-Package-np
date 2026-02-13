np.pairs <- function(y_vars, y_dat, ...) {
  if (missing(y_vars) || missing(y_dat))
    stop("'y_vars' and 'y_dat' are required")
  if (!is.data.frame(y_dat))
    y_dat <- as.data.frame(y_dat)
  if (is.null(names(y_vars)))
    names(y_vars) <- y_vars
  if (any(!y_vars %in% names(y_dat)))
    stop("all elements of 'y_vars' must be column names in 'y_dat'")

  pair_names <- expand.grid(y_vars, y_vars, stringsAsFactors = FALSE)
  pair_kerns <- lapply(seq_len(nrow(pair_names)), function(i) {
    y1 <- pair_names[i, 1]
    y2 <- pair_names[i, 2]
    if (y1 == y2) {
      npudens(tdat = y_dat[, y1], ...)
    } else {
      npreg(tydat = y_dat[, y2], txdat = y_dat[, y1], residuals = TRUE, ...)
    }
  })
  list(y_vars = y_vars, pair_names = pair_names, pair_kerns = pair_kerns)
}

np.pairs.plot <- function(pair_list) {
  if (length(pair_list) < 3)
    stop("pair_list must be created by np.pairs")

  pair_names <- pair_list[["pair_names"]]
  pair_kerns <- pair_list[["pair_kerns"]]
  y_vars <- pair_list[["y_vars"]]
  y_labels <- names(y_vars)

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  par(mfrow = c(length(y_vars), length(y_vars)), mar = c(4, 4, 2, 0))
  for (i in seq_len(nrow(pair_names))) {
    y1 <- pair_names[i, 1]
    y2 <- pair_names[i, 2]
    y1_label <- y_labels[match(y1, y_vars)]
    y2_label <- y_labels[match(y2, y_vars)]
    pair_i <- pair_kerns[[i]]
    if (y1 == y2) {
      ydat <- pair_i$eval[, 1]
      y_pred <- pair_i$dens
      y_ord <- order(ydat)
      y_pred <- y_pred[y_ord]
      ydat <- ydat[y_ord]
      plot.new()
      plot.window(xlim = range(ydat, na.rm = TRUE),
                  ylim = range(y_pred, na.rm = TRUE))
      polygon(c(ydat, rev(ydat)),
              c(y_pred, y_pred * 0),
              col = grey(0, 0.25), border = NA)
      axis(1)
      axis(2)
      mtext("Density Estimate", 2, line = 2.5, cex = 1)
      mtext(y1_label, 1, line = 2.5, cex = 1)
    } else {
      xdat <- pair_i$eval[, 1]
      y_pred <- pair_i$mean
      ydat <- y_pred + pair_i$resid
      x_ord <- order(xdat)
      y_pred <- y_pred[x_ord]
      ydat <- ydat[x_ord]
      xdat <- xdat[x_ord]
      plot.new()
      plot.window(xlim = range(xdat, na.rm = TRUE),
                  ylim = range(ydat, na.rm = TRUE))
      points(xdat, ydat, col = grey(0, 0.5), pch = 16, cex = 1.5)
      lines(xdat, y_pred, lwd = 2, col = rgb(0, 0, 1, 1))
      axis(1)
      axis(2)
      mtext(y2_label, 2, line = 2.5, cex = 1)
      mtext(y1_label, 1, line = 2.5, cex = 1)
    }
  }
  invisible(NULL)
}
