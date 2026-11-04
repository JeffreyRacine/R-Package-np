# Resolution is diagnostic: finite requested empirical bands remain available.
.np_band_resolution <- function(B, groups, alpha, active, bands, signed=FALSE,
                                 simultaneous.index=NULL) {
  index <- split(seq_along(groups),groups)
  rows <- list()
  for(g in seq_along(index)) {
    ii <- index[[g]];m <- length(ii)
    for(a in seq_along(alpha)) for(band in bands) {
      tail <- B*alpha[a]/if(band=="bonferroni")m else 1
      if(signed) tail <- tail/2
      order <- if(signed)1+floor(tail) else ceiling(B-tail)
      if(band=="simultaneous" && signed) order <- simultaneous.index[g,a]
      limited <- if(signed && band=="simultaneous")order==1 else tail<1
      rows[[length(rows)+1L]] <- data.frame(family=names(index)[g],band=band,
        alpha=alpha[a],B=B,m=m,nominal_tail_draws=tail,order_index=order,
        computed=TRUE,resolution_limited=any(active[ii]) && limited,
        calibration=if(signed)"rank" else "absolute",stringsAsFactors=FALSE)
    }
  }
  do.call(rbind,rows)
}

# Calibrated percentile families are active at every requested row. A row may
# use a different companion index at each point; report the weakest tail.
# nominal_tail_draws deliberately retains its nominal meaning.
.np_band_resolution_use_indices <- function(diagnostics, indices) {
  if (!is.data.frame(diagnostics) || !is.list(indices) ||
      length(indices) != nrow(diagnostics) ||
      !all(vapply(indices, function(x) is.numeric(x) && length(x) > 0L &&
        all(is.finite(x)) && all(x == floor(x)), logical(1))))
    stop("invalid calibrated percentile resolution indices", call. = FALSE)
  lower <- vapply(indices, min, numeric(1))
  upper <- vapply(indices, max, numeric(1))
  if (any(lower < 1 | upper > ceiling(diagnostics$B / 2)))
    stop("calibrated percentile resolution index outside draw bank", call. = FALSE)
  diagnostics$order_index <- unname(lower)
  diagnostics$resolution_limited <- unname(lower == 1)
  diagnostics$calibration <- "smoothed-percentile"
  diagnostics
}

.np_band_resolution_warn <- function(diagnostics) {
  if(is.list(diagnostics) && !is.data.frame(diagnostics)) {
    if(!length(diagnostics))return(invisible(NULL))
    diagnostics <- do.call(rbind,diagnostics)
  }
  if(is.null(diagnostics) || !any(diagnostics$resolution_limited))return(invisible(NULL))
  affected <- diagnostics[diagnostics$resolution_limited,,drop=FALSE]
  warning(structure(list(message=paste0("Empirical band quantiles have limited draw resolution (",
    paste(unique(affected$band),collapse=", "),
    "); all requested bands were computed. See the resolution diagnostics in the returned plot-data inference records."),call=NULL,
    diagnostics=affected),class=c("np_band_resolution_warning","warning","condition")))
  invisible(NULL)
}
