npksum <-
  function(...){
    mc <- match.call(expand.dots = FALSE)
    target <- .np_bw_dispatch_target(dots = mc$..., eval_env = parent.frame())
    UseMethod("npksum", target)
  }

.npRmpi_npksum_should_localize <- function(bws, dots = list()) {
  if (inherits(bws, "kbandwidth"))
    return(identical(bws$type, "adaptive_nn"))

  if (is.list(bws) && !is.null(bws$type))
    return(identical(bws$type, "adaptive_nn"))

  if (is.list(dots) && !is.null(dots$bwtype))
    return(identical(as.character(dots$bwtype)[1L], "adaptive_nn"))

  FALSE
}

npksum.formula <-
  function(formula, data, newdata, subset, na.action, ...){
    .npRmpi_require_active_slave_pool(where = "npksum()")
    if (.npRmpi_npksum_should_localize(NULL, list(...)) &&
        !isTRUE(getOption("npRmpi.local.regression.mode", FALSE)))
      return(.npRmpi_with_local_regression(.npRmpi_eval_without_dispatch(match.call(), parent.frame())))
    if (.npRmpi_autodispatch_active())
      return(.npRmpi_autodispatch_call(match.call(), parent.frame()))

    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(mf), nomatch = 0)
    mf <- mf[c(1,m)]
    mf[[1]] <- as.name("model.frame")
    mf.args <- as.list(mf[-1L])
    mf <- do.call(stats::model.frame, mf.args, envir = parent.frame())
    
    tydat <- model.response(mf)
    txdat <- mf[, attr(attr(mf, "terms"),"term.labels"), drop = FALSE]

    miss.new <- missing(newdata)
    if (!miss.new){
      tt <- delete.response(attr(mf,"terms"))
      umf.args <- list(formula = tt, data = newdata)
      umf <- do.call(stats::model.frame, umf.args, envir = parent.frame())
      exdat <- umf[, attr(attr(umf, "terms"),"term.labels"), drop = FALSE]
    }

    call_args <- list(txdat = txdat)
    if (!is.null(tydat))
      call_args$tydat <- tydat
    if (!miss.new)
      call_args$exdat <- exdat

    tobj <- do.call(npksum, c(call_args, list(...)))

    tobj$formula <- formula
    tobj$na.action <- attr(mf, "na.action")
    tobj$terms <- attr(mf,"terms")
    tobj
  }

npksum.integer <-
  function(...) { npksum.numeric(...) }

npksum.numeric <-
  function(bws,
           txdat = stop("training data 'txdat' missing"),
           tydat,
           exdat,
           weights,
           bandwidth.divide, compute.ocg, compute.score, kernel.pow,
           leave.one.out, operator, permutation.operator, return.kernel.weights,
           ...){
    dots <- list(...)
    return.derivative.kernel.weights <- isTRUE(dots$return.derivative.kernel.weights)
    dots$return.derivative.kernel.weights <- NULL
    .npRmpi_require_active_slave_pool(where = "npksum()")
    if (.npRmpi_npksum_should_localize(bws, list(...)) &&
        !isTRUE(getOption("npRmpi.local.regression.mode", FALSE)))
      return(.npRmpi_with_local_regression(.npRmpi_eval_without_dispatch(match.call(), parent.frame())))
    if (.npRmpi_autodispatch_active())
      return(.npRmpi_autodispatch_call(match.call(), parent.frame()))

    txdat <- toFrame(txdat)
    if (!missing(exdat)) {
      exdat_was_df <- is.data.frame(exdat)
      exdat_was_matrix <- is.matrix(exdat)
      exdat <- toFrame(exdat)
      # Preserve historical naming: `npksum.numeric()` used to forward `exdat=exdat`
      # into `npksum.default()`, so one-dimensional evaluation data typically carried
      # the column name "exdat" regardless of the caller's symbol.
      if (!exdat_was_df && !exdat_was_matrix && ncol(exdat) == 1L)
        names(exdat) <- "exdat"
    }

    kbw_args <- list(
      bw = bws,
      xdati = untangle(txdat),
      xnames = names(txdat)
    )
    if (!missing(tydat))
      kbw_args$ydati <- untangle(as.data.frame(tydat))

    kbw_args <- c(kbw_args, dots)
    tbw <- do.call(kbandwidth, kbw_args)

    call_args <- list(txdat = txdat, bws = tbw)
    if (!missing(tydat))
      call_args$tydat <- tydat
    if (!missing(exdat))
      call_args$exdat <- exdat
    if (!missing(weights))
      call_args$weights <- weights
    if (!missing(leave.one.out))
      call_args$leave.one.out <- leave.one.out
    if (!missing(kernel.pow))
      call_args$kernel.pow <- kernel.pow
    if (!missing(bandwidth.divide))
      call_args$bandwidth.divide <- bandwidth.divide
    if (!missing(operator))
      call_args$operator <- operator
    if (!missing(permutation.operator))
      call_args$permutation.operator <- permutation.operator
    if (!missing(compute.score))
      call_args$compute.score <- compute.score
    if (!missing(compute.ocg))
      call_args$compute.ocg <- compute.ocg
    if (!missing(return.kernel.weights))
      call_args$return.kernel.weights <- return.kernel.weights
    if (return.derivative.kernel.weights)
      call_args$return.derivative.kernel.weights <- return.derivative.kernel.weights

    do.call(npksum.default, call_args)
  }
    
npksum.default <- 
  function(bws,
           txdat = stop("training data 'txdat' missing"),
           tydat = NULL,
           exdat = NULL,
           weights = NULL,
           bandwidth.divide = FALSE,
           compute.ocg = FALSE,
           compute.score = FALSE,
           kernel.pow = 1.0,
           leave.one.out = FALSE,
           operator = names(ALL_OPERATORS),
           permutation.operator = names(PERMUTATION_OPERATORS),
           return.kernel.weights = FALSE,
           ...){
    dots <- list(...)
    return.derivative.kernel.weights <- isTRUE(dots$return.derivative.kernel.weights)
    .npRmpi_require_active_slave_pool(where = "npksum()")
    if (.npRmpi_npksum_should_localize(bws, list(...)) &&
        !isTRUE(getOption("npRmpi.local.regression.mode", FALSE)))
      return(.npRmpi_with_local_regression(.npRmpi_eval_without_dispatch(match.call(), parent.frame())))
    if (.npRmpi_autodispatch_active())
      return(.npRmpi_autodispatch_call(match.call(), parent.frame()))

    miss.ty <- missing(tydat)
    miss.ex <- missing(exdat)
    miss.weights <- missing(weights)

    uo.operators <- c("normal","convolution","integral")

    if(missing(operator))
      operator <- match.arg(operator, choices = names(ALL_OPERATORS))
    else
      operator <- match.arg(operator, choices = names(ALL_OPERATORS), several.ok = TRUE)

    permutation.operator <- match.arg(permutation.operator, choices = names(PERMUTATION_OPERATORS) )
    
    txdat = toFrame(txdat)

    if (!miss.ty && is.vector(tydat))
      tydat = as.matrix(tydat)

#    if(class(bws) != "kbandwidth")
    if(!isa(bws,"kbandwidth"))
        bws = kbandwidth(bws)

    if (!miss.ex){
      exdat = toFrame(exdat)

      if (! txdat %~% exdat )
        stop("'txdat' and 'exdat' are not similar data frames!")
    }

    if (length(bws$bw) != length(txdat))
      stop("length of bandwidth vector does not match number of columns of 'txdat'")

    if(length(operator) == 1)
      operator = rep(operator, length(txdat))
    
    if(length(operator) != length(txdat))
      stop("operator not specified for all variables")

    if(compute.score && compute.ocg)
      stop("compute.score and compute.ocg are mutually exclusive, and cannot be enabled simultaneously")
    if(!all(operator[bws$iuno | bws$iord] %in% uo.operators) && !compute.score)
      stop("unordered and ordered variables may only make use of 'normal', 'convolution' and 'integral' operator types")
    
    operator.num <- ALL_OPERATORS[operator]
    poperator.num <- PERMUTATION_OPERATORS[permutation.operator]
    
    if ((any(bws$icon) &&
         !all(vapply(txdat[, bws$icon, drop = FALSE], inherits, logical(1), c("integer", "numeric")))) ||
        (any(bws$iord) &&
         !all(vapply(txdat[, bws$iord, drop = FALSE], inherits, logical(1), "ordered"))) ||
        (any(bws$iuno) &&
         !all(vapply(txdat[, bws$iuno, drop = FALSE], inherits, logical(1), "factor"))))
      stop("supplied bandwidths do not match 'txdat' in type")

    if (!miss.ty && (nrow(txdat) != nrow(tydat)))
      stop("number of explanatory data 'txdat' and dependent data 'tydat' do not match")

    integer.pow = identical(as.double(as.integer(kernel.pow)),as.double(kernel.pow))
    
    if(!integer.pow)
      stop("'kernel.pow' is not an integer")

    if (!miss.ex && leave.one.out)
      stop("you may not specify 'leave.one.out = TRUE' and provide evaluation data")

    if (!miss.weights && !(is.matrix(weights) && nrow(weights) == nrow(txdat)))
      stop("improperly specified weight matrix")
    ## catch and destroy NA's

    keep.rows <- rep_len(TRUE, nrow(txdat))
    rows.omit <- attr(na.omit(txdat), "na.action")

    if (!miss.ty)
      rows.omit <- union(rows.omit, attr(na.omit(tydat), "na.action"))

    if (!miss.weights)
      rows.omit <- union(rows.omit, attr(na.omit(weights), "na.action"))

    
    if (length(rows.omit) > 0L)
      keep.rows[as.integer(rows.omit)] <- FALSE

    if (!any(keep.rows))
      stop("Data has no rows without NAs")

    txdat <- txdat[keep.rows,,drop = FALSE]

    if (!miss.ty)
      tydat <- tydat[keep.rows,, drop = FALSE]

    if (!miss.weights)
      weights <- weights[keep.rows,, drop = FALSE]

    if (!miss.ex){
      exdat = na.omit(exdat)
      
      if (nrow(exdat)==0)
        stop("Data has no rows without NAs")
    }

    ## end catch and destroy NA's

    tnrow = nrow(txdat)
    enrow = (if (miss.ex) tnrow else nrow(exdat))
    twncol = (if (miss.weights) 0 else ncol(weights))
    tyncol = (if (miss.ty) 0 else ncol(tydat))

    dim.in = c(twncol, tyncol, enrow)

    dim.out = c(twncol, tyncol, enrow)
    
    length.out = prod(dim.out[which(dim.out > 0)])

    has.permutation <- permutation.operator != "none"
    has.pksum <- has.permutation || compute.ocg

    if (has.pksum){
      npvar <- (if (has.permutation) bws$ncon else 0L) +
        (if (compute.score || compute.ocg) bws$nuno + bws$nord else 0L)
      p.length.out <- npvar*length.out
      p.dim.out <- c(dim.out, max(npvar, 0))      
    }
    else
      p.length.out <- 0

    
    ##dim.out = dim.out[dim.out > 1]

    ## re-assign levels in training and evaluation data to ensure correct
    ## conversion to numeric type.
    
    txdat <- adjustLevels(txdat, bws$xdati, allowNewCells = TRUE)
    
    if (!miss.ex)
      exdat <- adjustLevels(exdat, bws$xdati, allowNewCells = TRUE)

    if (!miss.ex)
      npKernelBoundsCheckEval(exdat, bws$icon, bws$ckerlb, bws$ckerub, argprefix = "cker")

    ## grab the evaluation data before it is converted to numeric
    if(miss.ex)
      teval <- txdat
    else
      teval <- exdat

    ## put the unordered, ordered, and continuous data in their own objects
    ## data that is not a factor is continuous.
    
    txdat = toMatrix(txdat)

    tuno = txdat[, bws$iuno, drop = FALSE]
    tcon = txdat[, bws$icon, drop = FALSE]
    tord = txdat[, bws$iord, drop = FALSE]

    if (!miss.ex){
      exdat = toMatrix(exdat)

      euno = exdat[, bws$iuno, drop = FALSE]
      econ = exdat[, bws$icon, drop = FALSE]
      eord = exdat[, bws$iord, drop = FALSE]

    } else {
      euno = data.frame()
      eord = data.frame()
      econ = data.frame()
    }

    nkw <- (if (return.kernel.weights) tnrow*enrow else 0)

    need.pkw <- isTRUE(return.derivative.kernel.weights) &&
      return.kernel.weights &&
      has.pksum &&
      (p.length.out > 0L)

    return.names <- c("ksum", "kernel.weights", "p.ksum")
    if (need.pkw)
      return.names <- c(return.names, "p.kernel.weights")
      
	    myopti = list(
      num_obs_train = tnrow,
      num_obs_eval = enrow,
      num_uno = bws$nuno,
      num_ord = bws$nord,
      num_con = bws$ncon,
      int_LARGE_SF = SF_ARB,
      BANDWIDTH_reg_extern = switch(bws$type,
        fixed = BW_FIXED,
        generalized_nn = BW_GEN_NN,
        adaptive_nn = BW_ADAP_NN),
      int_MINIMIZE_IO=if (isTRUE(getOption("np.messages"))) IO_MIN_FALSE else IO_MIN_TRUE, 
      kerneval = switch(bws$ckertype,
        gaussian = CKER_GAUSS + bws$ckerorder/2 - 1,
        epanechnikov = CKER_EPAN + bws$ckerorder/2 - 1,
        uniform = CKER_UNI,
        "truncated gaussian" = CKER_TGAUSS),
      ukerneval = switch(bws$ukertype,
        aitchisonaitken = UKER_AIT,
        liracine = UKER_LR),
      okerneval = switch(bws$okertype,
        wangvanryzin = OKER_WANG,
        liracine = OKER_LR,
        nliracine = OKER_NLR,
        "racineliyan" = OKER_RLY),
      miss.ex = miss.ex,
      leave.one.out = leave.one.out, 
      bandwidth.divide = bandwidth.divide,
      mcv.numRow = attr(bws$xmcv, "num.row"),
      wncol = dim.in[1],
      yncol = dim.in[2],
      int_do_tree = if (isTRUE(getOption("np.tree"))) DO_TREE_YES else DO_TREE_NO,
      return.kernel.weights = return.kernel.weights,
      permutation.operator = poperator.num,
      compute.score = compute.score,
      compute.ocg = compute.ocg,
      suppress.parallel = isTRUE(getOption("npRmpi.local.regression.mode", FALSE)))

	    cker.bounds.c <- npKernelBoundsMarshal(bws$ckerlb[bws$icon], bws$ckerub[bws$icon])
    
   asDouble <- function(data){
	   if (is.null(data)){
	 	   result <- as.double(0.0)
	   }
	   else {
		   result <- as.double(data)
	   }
	   return(result)
   }

    myout <-
      .Call("C_np_kernelsum",
            asDouble(tuno), asDouble(tord), asDouble(tcon),
            asDouble(tydat), asDouble(weights),
            asDouble(euno), asDouble(eord), asDouble(econ),
            as.double(c(bws$bw[bws$icon], bws$bw[bws$iuno], bws$bw[bws$iord])),
            as.double(bws$xmcv), as.double(attr(bws$xmcv, "pad.num")),
            as.integer(c(operator.num[bws$icon], operator.num[bws$iuno], operator.num[bws$iord])),
            as.integer(myopti), as.double(kernel.pow),
            as.integer(length.out),
            as.integer(p.length.out),
            as.integer(nkw),
            as.double(cker.bounds.c$lb),
            as.double(cker.bounds.c$ub),
            PACKAGE="npRmpi")[return.names]

    if (dim.out[1] > 1){
      dim(myout[["ksum"]]) <- dim.out
      if (dim.out[2] < 2)
        dim(myout[["ksum"]]) <- dim(myout[["ksum"]])[-2]
    } else if (max(dim.out) > 1) {
      dim(myout[["ksum"]]) <- dim.out[dim.out > 1]
      if (miss.weights)
        myout[["ksum"]] <- aperm(myout[["ksum"]])
    }

    if(return.kernel.weights){
      kw <- matrix(data = myout[["kernel.weights"]], nrow = tnrow, ncol = enrow)
    } else {
      kw <- NULL
    }

    if(need.pkw){
      raw.pkw <- myout[["p.kernel.weights"]]
      if(length(raw.pkw) > 0L){
        if(npvar == 1L){
          p.kw <- matrix(data = raw.pkw, nrow = tnrow, ncol = enrow)
        } else {
          p.kw <- array(data = raw.pkw, dim = c(tnrow, enrow, npvar))
        }
      } else {
        p.kw <- NULL
      }
    } else {
      p.kw <- NULL
    }

    if(has.pksum && (p.length.out > 0)) {
      dim.p <- p.dim.out[which(p.dim.out > 1)]
      if(length(dim.p) == 0) dim.p <- 1

      p.myout <- matrix(data = myout[["p.ksum"]], ncol = npvar)

      ip <- integer(0)

      if(has.permutation && (bws$ncon > 0))
        ip <- which(bws$icon)

      if(compute.ocg || compute.score){
        if(bws$nuno > 0)
          ip <- c(ip,which(bws$iuno))

        if(bws$nord > 0)
          ip <- c(ip,which(bws$iord))
      }

      p.myout[,ip] <- p.myout

      p.myout <- array(data = as.vector(p.myout), dim = dim.p)
    } else {
      p.myout <- NULL
    }
    return( npkernelsum(bws = bws,
                        eval = teval,
                        ksum = myout[["ksum"]],
                        kw = kw,
                        p.ksum = p.myout,
                        p.kw = p.kw,
                        ntrain = tnrow, trainiseval = miss.ex) )

  }
