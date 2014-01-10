npksum <-
  function(...){
    args = list(...)
    if (is(args[[1]],"formula"))
      UseMethod("npksum",args[[1]])
    else if (!is.null(args$formula))
      UseMethod("npksum",args$formula)
    else
      UseMethod("npksum",args[[which(names(args)=="bws")[1]]])
  }

npksum.formula <-
  function(formula, data, newdata, subset, na.action, ...){

    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(mf), nomatch = 0)
    mf <- mf[c(1,m)]
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    
    tydat <- model.response(mf)
    txdat <- mf[, attr(attr(mf, "terms"),"term.labels"), drop = FALSE]

    if (!(miss.new <- missing(newdata))){
      tt <- delete.response(attr(mf,"terms"))
      umf <- emf <- model.frame(tt, data = newdata)
      exdat <- emf[, attr(attr(emf, "terms"),"term.labels"), drop = FALSE]
    }

    tobj <- eval(parse(text=paste("npksum(txdat = txdat,",
                         ifelse(!is.null(tydat), "tydat = tydat,", ""),
                         ifelse(!miss.new,"exdat = exdat,",""), "...)")))

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
           leave.one.out, kernel.pow, bandwidth.divide,
           operator, permutation.operator, compute.score, compute.ocg, return.kernel.weights,
           ...){

    txdat <- toFrame(txdat)

    
    tbw <- eval(parse(text=paste("kbandwidth(bw = bws,",
                        "xdati = untangle(txdat),",
                        ifelse(missing(tydat),"","ydati = untangle(as.data.frame(tydat)),"),
                        "xnames = names(txdat), ...)")))

    mc.names <- names(match.call(expand.dots = FALSE))
    margs <- c("tydat", "exdat", "weights", "leave.one.out", "kernel.pow", "bandwidth.divide",
               "operator", "permutation.operator", "compute.score", "compute.ocg", "return.kernel.weights")
    m <- match(margs, mc.names, nomatch = 0)
    any.m <- any(m != 0)


    eval(parse(text=paste("npksum.default(txdat=txdat, bws=tbw",
                          ifelse(any.m, ",",""),
                          paste(mc.names[m], ifelse(any.m,"=",""), mc.names[m], collapse=", "),
                          ")")))
  }
    
npksum.default <- 
  function(bws,
           txdat = stop("training data 'txdat' missing"),
           tydat = NULL,
           exdat = NULL,
           weights = NULL,
           leave.one.out = FALSE,
           kernel.pow = 1.0,
           bandwidth.divide = FALSE,
           operator = names(ALL_OPERATORS),
           permutation.operator = names(PERMUTATION_OPERATORS),
           compute.score = FALSE,
           compute.ocg = FALSE,
           return.kernel.weights = FALSE,
           ...){

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

    if(class(bws) != "kbandwidth")
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
    
    ccon = unlist(lapply(txdat[,bws$icon,drop=FALSE],class))
    if ((any(bws$icon) && !all((ccon == class(integer(0))) | (ccon == class(numeric(0))))) ||
        (any(bws$iord) && !all(unlist(lapply(txdat[,bws$iord, drop=FALSE],class)) ==
                               class(ordered(0)))) ||
        (any(bws$iuno) && !all(unlist(lapply(txdat[,bws$iuno, drop=FALSE],class)) ==
                               class(factor(0)))))
      stop("supplied bandwidths do not match 'txdat' in type")

    if (!miss.ty && (nrow(txdat) != nrow(tydat)))
      stop("number of explanatory data 'txdat' and dependent data 'tydat' do not match")

    integer.pow = identical(as.double(as.integer(kernel.pow)),as.double(kernel.pow))
    
    if(!integer.pow)
      stop("'kernel.pow' is not an integer")

    if (!miss.ex & leave.one.out)
      stop("you may not specify 'leave.one.out = TRUE' and provide evaluation data")

    if (!miss.weights && !(is.matrix(weights) && nrow(weights) == nrow(txdat)))
      stop("improperly specified weight matrix")
    ## catch and destroy NA's

    goodrows = 1:dim(txdat)[1]
    rows.omit = attr(na.omit(txdat), "na.action")

    if (!miss.ty)
      rows.omit = union(rows.omit, attr(na.omit(tydat), "na.action"))

    if (!miss.weights)
      rows.omit = union(rows.omit, attr(na.omit(weights), "na.action"))

    
    goodrows[rows.omit] = 0

    if (all(goodrows==0))
      stop("Data has no rows without NAs")

    txdat = txdat[goodrows,,drop = FALSE]

    if (!miss.ty)
      tydat = tydat[goodrows,, drop = FALSE]

    if (!miss.weights)
      weights = weights[goodrows,, drop = FALSE]

    if (!miss.ex){
      exdat = na.omit(exdat)
      
      if (nrow(exdat)==0)
        stop("Data has no rows without NAs")
    }

    ## end catch and destroy NA's

    tnrow = nrow(txdat)
    enrow = ifelse(miss.ex, tnrow, nrow(exdat))
    twncol = ifelse(miss.weights, 0, ncol(weights))
    tyncol = ifelse(miss.ty, 0, ncol(tydat))

    dim.in = c(twncol, tyncol, enrow)

    dim.out = c(max(ncol(weights),0), max(ncol(tydat),0), enrow)
    
    length.out = prod(dim.out[which(dim.out > 0)])

    if((permutation.operator != "none") || compute.ocg){
      npvar <- ifelse((permutation.operator != "none"), bws$ncon, 0) + ifelse(compute.score | compute.ocg, bws$nuno + bws$nord, 0)
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

    nkw <- ifelse(return.kernel.weights, tnrow*enrow, 0)

    return.names <- c("ksum","kernel.weights","p.ksum")
      
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
      int_MINIMIZE_IO=ifelse(options('np.messages'), IO_MIN_FALSE, IO_MIN_TRUE), 
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
        nliracine = OKER_NLR),
      miss.ex = miss.ex,
      leave.one.out = leave.one.out, 
      bandwidth.divide = bandwidth.divide,
      mcv.numRow = attr(bws$xmcv, "num.row"),
      wncol = dim.in[1],
      yncol = dim.in[2],
      int_do_tree = ifelse(options('np.tree'), DO_TREE_YES, DO_TREE_NO),
      return.kernel.weights = return.kernel.weights,
      permutation.operator = poperator.num,
      compute.score = compute.score,
      compute.ocg = compute.ocg)
    


    myout <- 
      .C("np_kernelsum",
         as.double(tuno), as.double(tord), as.double(tcon),
         as.double(tydat), as.double(weights),
         as.double(euno),  as.double(eord),  as.double(econ), 
         as.double(c(bws$bw[bws$icon],bws$bw[bws$iuno],bws$bw[bws$iord])),
         as.double(bws$xmcv), as.double(attr(bws$xmcv, "pad.num")),
         as.integer(c(operator.num[bws$icon],operator.num[bws$iuno],operator.num[bws$iord])),
         as.integer(myopti), as.double(kernel.pow),
         ksum = double(length.out),
         p.ksum = double(p.length.out),
         kernel.weights = double(nkw),
         PACKAGE="np" )[return.names]

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

    if(((permutation.operator != "none") || compute.ocg) && (p.length.out > 0)) {
      dim.p <- p.dim.out[which(p.dim.out > 1)]
      if(length(dim.p) == 0) dim.p <- 1

      p.myout <- matrix(data = myout[["p.ksum"]], ncol = npvar)

      ip <- integer(0)

      if((permutation.operator != "none") && (bws$ncon > 0))
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
                        ntrain = tnrow, trainiseval = miss.ex) )

  }
