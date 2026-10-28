# Acquisition policy only; multivariate variance/studentization is unchanged.
.np_bias_main_capture_supported <- function(state, law) {
  if (!is.character(law) || length(law)!=1L ||
      !law %in% c("gaussian","mammen") || state$d<2L ||
      length(state$factor.cols) || isTRUE(state$pooled) ||
      isTRUE(state$categorical.only) || !is.null(state$adaptation) ||
      !is.null(state$mixed) ||
      !.np_plot_regression_capture_supported(state$bws) ||
      (!is.null(state$gradient.order) && any(state$gradient.order!=1L)))
    return(FALSE)
  reg <- npValidatedConditionalRegSpec(state$bws,
    where="captured main bias rows",ncon.field="ncon")
  identical(reg[["basis.engine"]],"glp")
}

# Private, explicit integration seam for an already captured original fit.
# Acquisition belongs to the existing plot evaluator. Replay never refits.
.np_bias_capture_prepare <- function(spec, fit, grid, targets = "mean") {
  if (isTRUE(spec$pooled) || isTRUE(spec$categorical.only) ||
      length(spec$factor.cols) || !is.null(spec$adaptation) ||
      !is.null(spec$mixed) || !.np_plot_regression_capture_supported(spec$bws))
    stop("captured local bias prototype requires an inactive continuous fixed Gaussian2 request",call.=FALSE)
  request <- fit[["fit.request",exact=TRUE]]
  if (!is.list(request) || !identical(request$exdat,grid) ||
      !is.numeric(request$ydat) || !identical(as.double(request$ydat),spec$y))
    stop("captured local bias request does not match its grid or response",call.=FALSE)
  x <- request$xdat
  if (!is.data.frame(x) || !identical(names(x),names(spec$x)) ||
      nrow(x)!=nrow(spec$x))
    stop("captured local bias request does not match its training data",call.=FALSE)
  for (j in seq_len(spec$d))
    x[[spec$numeric.cols[j]]] <- (x[[spec$numeric.cols[j]]]-spec$offset[j])/spec$scale[j]
  h <- unname(request$bws[["bandwidth"]][["x"]])
  h[spec$numeric.cols] <- h[spec$numeric.cols]/spec$scale
  expected <- .np_bias_regbw(request$bws,x,h)
  if (!identical(x,spec$x) ||
      !.np_plot_regression_capture_same_bws(expected,spec$bws))
    stop("captured local bias request does not match the fixed model",call.=FALSE)
  valid <- c("mean",paste0("d",seq_len(spec$d)))
  if (!is.character(targets) || !length(targets) || anyNA(targets) ||
      anyDuplicated(targets) || any(!targets %in% valid))
    stop("invalid captured local bias target family",call.=FALSE)
  axes <- setNames(vapply(targets,function(t)
    if(t=="mean")0L else as.integer(substring(t,2L)),integer(1)),targets)
  if (any(vapply(axes,function(j)j>0L && .np_bias_derivative_order(spec,j)!=1L,logical(1))))
    stop("captured local bias prototype admits first derivatives only",call.=FALSE)
  operators <- setNames(lapply(axes,function(j) {
    nu <- NULL
    if(j) {nu <- integer(spec$d);nu[j] <- 1L}
    .np_plot_regression_capture_operator(fit,request$bws,request$xdat,
      request$ydat,grid,nu,NULL,64*1024^2)$operator
  }),targets)
  ridge <- fit[["fit.inference",exact=TRUE]][[3L]]
  reg <- npValidatedConditionalRegSpec(spec$bws,where="bias main operator",ncon.field="ncon")
  # The complete request was validated above. Bind admission to this exact
  # correction state without selecting a variance or calibration policy.
  certificate <- if (!is.null(spec$local.inference)) TRUE else
    structure(list(state=spec), class="np_bias_main_capture_certificate")
  # The immutable complete grid is indexed by occurrence, never value-matched.
  rows <- function(index,map) {
    index <- .np_plot_lsq_gate_indices(index,nrow(grid))
    if (length(index)>25L || !is.data.frame(map) || nrow(map)<1L ||
        nrow(map)>500L || !is.numeric(map$row) || anyNA(map$row) ||
        any(map$row<1L | map$row>length(index) | map$row!=floor(map$row)) ||
        anyNA(map$target) || any(!map$target %in% targets))
      stop("invalid bounded captured local bias tile",call.=FALSE)
    L <- matrix(0,nrow(map),nrow(spec$x));diagnostic <- list()
    for (target in unique(map$target)) {
      rr <- which(map$target==target);axis <- axes[[target]]
      nu <- NULL
      if(axis) {nu <- integer(spec$d);nu[axis] <- 1L}
      for (first in seq.int(1L,length(rr),by=16L)) {
        jj <- rr[seq.int(first,min(length(rr),first+15L))]
        ii <- index[map$row[jj]]
        H <- operators[[target]](grid[ii,,drop=FALSE],nu,ii)
        H <- .np_bias_main_admit(spec,H,ridge[ii],reg,
          captured=certificate)
        # Auxiliary correction operates on standardized coordinates. Convert
        # here; the existing final conversion restores physical units once.
        L[jj,] <- H * if(axis)spec$scale[axis] else 1
        diagnostic[[length(diagnostic)+1L]] <- data.frame(row=map$row[jj],
          target=target,frame=1L,ridge=attr(H,"ridge.used",exact=TRUE),
          stringsAsFactors=FALSE)
      }
    }
    attr(L,"main.regularization") <- do.call(rbind,diagnostic)
    .np_bias_regression_rows(spec,grid[index,,drop=FALSE],map,original=L,
      captured=certificate)
  }
  structure(list(state=spec,grid=grid,targets=targets,rows=rows),
    class="np_bias_capture_prototype")
}

.np_bias_capture_validate <- function(capture, state, grid, targets) {
  if (!inherits(capture,"np_bias_capture_prototype") ||
      !identical(capture$state,state) || !identical(capture$grid,grid) ||
      !identical(capture$targets,targets) || !is.function(capture$rows))
    stop("captured local bias state does not match inference",call.=FALSE)
  invisible(NULL)
}

