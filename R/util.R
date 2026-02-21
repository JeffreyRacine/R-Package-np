## This function will compute the cumulative integral at each sample
## realization using the trapezoidal rule and the cumsum function as
## we need to compute this in a computationally efficient manner.

# integrate.trapezoidal <- function(x,y) {
#   n <- length(x)
#   rank.x <- rank(x)
#   order.x <- order(x)
#   y <- y[order.x]
#   x <- x[order.x]
#   int.vec <- numeric(length(x))
#   ## Use a correction term at the boundary: -cx^2/12*(f'(b)-f'(a)),
#   ## check for NaN case
#   cx  <- x[2]-x[1]
#   ca <- (y[2]-y[1])/cx
#   cb <- (y[n]-y[n-1])/cx
#   cf <- cx^2/12*(cb-ca)
#   if(!is.finite(cf)) cf <- 0
#   int.vec[1] <- 0
#   int.vec[2:n] <- cumsum((x[2:n]-x[2:n-1])*(y[2:n]+y[2:n-1])/2)
#   return(int.vec[rank.x]-cf)
# }

# Benches 1.3-1.7 times faster (mean/median) produces identical results

integrate.trapezoidal <- function(x, y) {
  n <- length(x)
  order.x <- order(x)
  x <- x[order.x]
  y <- y[order.x]
  dx <- diff(x)
  dy <- diff(y)
  cx <- dx[1]
  ca <- dy[1] / cx
  cb <- dy[n - 1] / cx
  cf <- cx^2 / 12 * (cb - ca)
  if (!is.finite(cf)) cf <- 0
  int.vec <- c(0, cumsum(dx * (y[-n] + y[-1]) / 2))
  int.vec <- int.vec - cf
  int.vec[order(order.x)]  # inverse permutation == rank(x) when no ties
}

## No Zero Denominator, used in C code for kernel estimation...

## Original, better, best

# NZD <- function(a) {
#   ifelse(a<0,pmin(-.Machine$double.eps,a),pmax(.Machine$double.eps,a))
# }
  
# NZD <- function(a) {
#   if(length(a) == 1) {
#     if(is.na(a)) return(a)
#     if(a < 0) return(min(-.Machine$double.eps, a))
#     return(max(.Machine$double.eps, a))
#   }
#   ifelse(a<0,pmin(-.Machine$double.eps,a),pmax(.Machine$double.eps,a))
# }

# Benches 5.5-7.9 times faster (mean/median) produces identical results

NZD <- function(a) {
  eps <- .Machine$double.eps
  if (length(a) == 1) {
    if (a >= 0) {
      if (a < eps) return(eps)
    } else {
      if (a > -eps) return(-eps)
    }
    return(a)
  }
  idx <- which(abs(a) < eps)
  if (length(idx) > 0)
    a[idx] <- ifelse(a[idx] >= 0, eps, -eps)
  a
}

## New function when only positive values are guaranteed, better, best

# NZD_pos <- function(a) {
#   if(length(a) == 1) {
#     if(is.na(a)) return(a)
#     return(max(.Machine$double.eps, a))
#   }
#   pmax(.Machine$double.eps,a)
# }

# Benches 1.3-1.8 times faster (mean/median) produces identical results

NZD_pos <- function(a) {
  eps <- .Machine$double.eps
  if (length(a) == 1)
    return(if (a < eps) eps else a)
  idx <- which(a < eps)
  if (length(idx) > 0)
    a[idx] <- eps
  a
}

npValidateGlpDegree <- function(regtype, degree, ncon, argname = "degree") {
  degree.max <- 12L

  if (!identical(regtype, "lp"))
    return(NULL)

  if (is.null(degree))
    degree <- rep.int(1L, ncon)

  if (!length(degree) && ncon == 0L)
    return(integer(0))

  if (length(degree) != ncon)
    stop(sprintf("%s must have one entry per continuous predictor (%d expected, got %d)",
                 argname, ncon, length(degree)))

  if (any(!is.finite(degree)))
    stop(sprintf("%s must contain finite non-negative integers", argname))

  if (any(degree < 0))
    stop(sprintf("%s must contain finite non-negative integers", argname))

  if (any(degree != floor(degree)))
    stop(sprintf("%s must contain finite non-negative integers", argname))

  if (any(degree > degree.max))
    stop(sprintf("%s must contain finite non-negative integers in [0,%d]",
                 argname, degree.max))

  as.integer(degree)
}

npValidateGlpBernstein <- function(regtype, bernstein.basis, argname = "bernstein.basis") {
  if (!identical(regtype, "lp"))
    return(FALSE)

  if (is.null(bernstein.basis))
    bernstein.basis <- FALSE

  if (!is.logical(bernstein.basis) || length(bernstein.basis) != 1L || is.na(bernstein.basis))
    stop(sprintf("%s must be TRUE or FALSE", argname))

  isTRUE(bernstein.basis)
}

npValidateLpBasis <- function(regtype, basis, argname = "basis") {
  if (!identical(regtype, "lp"))
    return("glp")

  if (is.null(basis))
    basis <- "glp"

  basis <- match.arg(as.character(basis), c("glp", "additive", "tensor"))
  basis
}

npLpBasisCode <- function(basis) {
  switch(tolower(ifelse(is.null(basis) || !length(basis), "glp", basis)),
         additive = 0L,
         glp = 1L,
         tensor = 2L,
         1L)
}

npValidateGlpGradientOrder <- function(regtype,
                                       gradient.order,
                                       ncon,
                                       argname = "gradient.order") {
  degree.max <- 12L

  if (!identical(regtype, "lp"))
    return(NULL)

  if (is.null(gradient.order))
    gradient.order <- rep.int(1L, ncon)

  if (!length(gradient.order) && ncon == 0L)
    return(integer(0))

  if (length(gradient.order) == 1L && ncon > 1L)
    gradient.order <- rep.int(gradient.order, ncon)

  if (length(gradient.order) != ncon)
    stop(sprintf("%s must have one entry per continuous predictor (%d expected, got %d)",
                 argname, ncon, length(gradient.order)))

  if (any(!is.finite(gradient.order)))
    stop(sprintf("%s must contain finite positive integers", argname))

  if (any(gradient.order <= 0))
    stop(sprintf("%s must contain finite positive integers", argname))

  if (any(gradient.order != floor(gradient.order)))
    stop(sprintf("%s must contain finite positive integers", argname))

  if (any(gradient.order > degree.max))
    stop(sprintf("%s must contain integers in [1,%d]", argname, degree.max))

  as.integer(gradient.order)
}

npCheckRegressionDesignCondition <- function(reg.code,
                                             xcon,
                                             basis = "glp",
                                             degree = NULL,
                                             bernstein.basis = FALSE,
                                             where = "npregbw") {
  kappa.warn <- 1e8
  kappa.stop <- 1e12

  if (!(reg.code %in% c(REGTYPE_LL, REGTYPE_GLP)))
    return(invisible(NULL))

  xcon <- as.data.frame(xcon)
  n <- nrow(xcon)
  if (is.null(n) || n <= 0L)
    return(invisible(NULL))

  B <- if (identical(reg.code, REGTYPE_GLP)) {
    if (is.null(degree))
      stop(sprintf("%s: LP degree vector missing for design-conditioning check", where))
    W.lp(xdat = xcon,
         degree = degree,
         basis = basis,
         bernstein.basis = isTRUE(bernstein.basis))
  } else {
    cbind(1, as.matrix(xcon))
  }

  p <- ncol(B)
  if (is.null(p) || p <= 0L)
    return(invisible(NULL))

  sv <- suppressWarnings(tryCatch(svd(B, nu = 0L, nv = 0L)$d, error = function(e) NULL))
  if (is.null(sv) || !length(sv))
    stop(sprintf("%s: unable to compute singular values for design-conditioning check", where))
  tol.rank <- max(dim(B)) * max(sv) * .Machine$double.eps
  r <- sum(sv > tol.rank)
  if (r < p) {
    stop(sprintf("%s: regression design matrix is rank deficient (rank=%d < p=%d). Reduce polynomial degree or remove collinear continuous predictors.",
                 where, r, p))
  }

  kB <- suppressWarnings(tryCatch(kappa(B), error = function(e) Inf))
  if (!is.finite(kB))
    kB <- Inf

  if (kB > kappa.stop) {
    stop(sprintf("%s: regression design matrix is severely ill-conditioned (kappa(B)=%.3e > %.1e). Reduce polynomial degree or remove collinear continuous predictors.",
                 where, kB, kappa.stop))
  }

  if (kB > kappa.warn) {
    warning(sprintf("%s: regression design matrix is ill-conditioned (kappa(B)=%.3e > %.1e). Estimation may rely heavily on ridging; consider lower degree or less collinear predictors.",
                    where, kB, kappa.warn),
            call. = FALSE, immediate. = TRUE)
  }

  invisible(NULL)
}

npRegtypeToC <- function(regtype, degree, ncon, context = "npreg") {
  if (identical(regtype, "lc"))
    return(list(code = REGTYPE_LC, degree = NULL))

  if (identical(regtype, "ll"))
    return(list(code = REGTYPE_LL, degree = NULL))

  degree <- npValidateGlpDegree(regtype, degree, ncon)

  if ((ncon == 0L) || all(degree == 0L))
    return(list(code = REGTYPE_LC, degree = degree))

  if (all(degree == 1L))
    return(list(code = REGTYPE_LL, degree = degree))

  list(code = REGTYPE_GLP, degree = degree)
}

npRejectLegacyLpArgs <- function(dotnames, where = "npreg") {
  if (is.null(dotnames) || !length(dotnames))
    return(invisible(NULL))
  bad <- intersect(dotnames, c("glp.degree", "glp.bernstein", "glp.basis"))
  if (length(bad))
    stop(sprintf("%s: legacy arguments %s are no longer supported; use basis, degree and bernstein.basis",
                 where, paste(sprintf("'%s'", bad), collapse = ", ")))
  invisible(NULL)
}

dim_basis <- function(basis = c("additive", "glp", "tensor"),
                      kernel = TRUE,
                      degree = NULL,
                      segments = NULL,
                      include = NULL,
                      categories = NULL) {
  basis <- match.arg(basis)

  if (is.null(degree))
    degree <- integer(0)
  if (is.null(segments))
    segments <- rep.int(1L, length(degree))

  degree <- as.integer(degree)
  segments <- as.integer(segments)

  if (length(degree) != length(segments))
    stop("degree and segments must have the same length")
  if (any(!is.finite(degree)) || any(degree < 0L))
    stop("degree must contain finite non-negative integers")
  if (any(!is.finite(segments)) || any(segments <= 0L))
    stop("segments must contain finite positive integers")

  if (is.null(include))
    include <- integer(0)
  if (is.null(categories))
    categories <- integer(0)

  include <- as.integer(include)
  categories <- as.integer(categories)

  if (length(include) != length(categories))
    stop("include and categories must have the same length")
  if (any(!is.finite(include)) || any(include < 0L))
    stop("include must contain finite non-negative integers")
  if (any(!is.finite(categories)) || any(categories < 0L))
    stop("categories must contain finite non-negative integers")

  basis.code <- switch(basis, additive = 0L, glp = 1L, tensor = 2L)

  .C("np_dim_basis",
     basis_code = as.integer(basis.code),
     kernel = as.integer(isTRUE(kernel)),
     degree = degree,
     segments = segments,
     k = as.integer(length(degree)),
     include = include,
     categories = categories,
     ninclude = as.integer(length(include)),
     result = double(1L),
     PACKAGE = "npRmpi")$result
}

## Function to test for monotone increasing vector

is.monotone.increasing <- function(x) {
  ## Sorted and last value > first value
  !is.unsorted(x) && x[length(x)] > x[1]
}
  
## what is a badord? ... an ordered factor of numeric values to treat
## them properly one must preserve the numeric value, ie. scale not
## just their sorted order Actually, the ord/badord paradigm must go,
## in place of levels caching

matrix.sd <- function(x, na.rm=FALSE) {
  if(is.matrix(x)) apply(x, 2, sd, na.rm=na.rm)
  else if(is.vector(x)) sd(x, na.rm=na.rm)
  else if(is.data.frame(x)) sapply(x, sd, na.rm=na.rm)
  else sd(as.vector(x), na.rm=na.rm)
}

npseed <- function(seed){
  .C("np_set_seed",as.integer(abs(seed)), PACKAGE = "npRmpi")
  invisible()
}

erf <- function(z) { 2 * pnorm(z*sqrt(2)) - 1 }

nptgauss <- function(b){

  rel.tol <- sqrt(.Machine$double.eps)

  b.max <- sqrt(-2*log(.Machine$double.eps))

  if((b < 0) || (b > b.max))
    stop(paste("b must be between 0 and",b.max))

  alpha <- 1.0/(pnorm(b)-pnorm(-b)-2*b*dnorm(b))

  tgauss <- function(z)
    ifelse(abs(z) >= b, 0.0, alpha*(dnorm(z) - dnorm(b)))

  c0 <- alpha*dnorm(b)

  k <- integrate(f = function(z) { tgauss(z)^2 }, -b, b)$value
  k2 <- integrate(f = function(z) { z^2*tgauss(z) }, -b, b)$value
  k22 <- integrate(f = function(z) { (z*tgauss(z))^2 }, -b, b)$value
  km <- integrate(f = function(z) { tgauss(z-0.5)*tgauss(z+0.5) }, -b+0.5, b-0.5)$value

  a0 <- (0.5 + 2*b*c0)/integrate(f = function(z){ erf(z/2 + b)*exp(-0.25*z^2) }, -2*b, 0)$value
  a2 <- (c0 + k - a0*erf(b))/erf(b/sqrt(2))
  a1 <- -(a2*erf(b/sqrt(2)) + c0)/(2*b)

  int.kernels[CKER_TGAUSS + 1] <- k
  
  invisible(.C("np_set_tgauss2",as.double(c(b, alpha, c0, a0, a1, a2, k, k2, k22, km)), PACKAGE = "npRmpi"))

}

numNotIn <- function(x){
  while(is.element(num <- rnorm(1),x)){}
  num
}

dlev <- function(x){
  if(is.ordered(x))
    x.dlev <- suppressWarnings(as.numeric(levels(x)))
  if (!is.ordered(x) || any(is.na(x.dlev)))
    x.dlev <- as.numeric(1:nlevels(x))
  x.dlev
}

isNum <- function(x){
  return(!any(is.na(suppressWarnings(as.numeric(x)))))
}

untangle <- function(frame){
  if (is.null(frame))
    return(NULL)
  
  iord <- unlist(lapply(frame, is.ordered))
  iuno <- unlist(lapply(frame, is.factor)) & !iord
  icon <- unlist(lapply(frame, is.numeric))

  if(!all(iord | iuno | icon)) 
    stop("non-allowed data type in data frame")

  inumord <-
    suppressWarnings(unlist(lapply(frame,
    function (z) {
      is.ordered(z) && is.numeric(tryCatch(as.numeric(levels(z)), warning = function (y) {
        FALSE }))
    })))

  all.lev <- lapply(frame, function(y){
    t.ret <- NULL
    if (is.factor(y))
      t.ret <- levels(y)
    t.ret
  })

  all.ulev <- lapply(frame, function (y) {
    t.ret <- NULL
    if (is.factor(y))
      t.ret <- sort(unique(y))
    t.ret
  })

  all.dlev <- lapply(frame, function (y) {
    t.ret <- NULL
    if (is.factor(y))
      t.ret <- dlev(y)
    t.ret
  })

  all.nlev <- lapply(frame, function(y) {
    t.ret <- NULL
    if (is.factor(y))
      t.ret <- nlevels(y)
    t.ret
  })

  all.min <- lapply(frame, function(y) {
    t.ret <- NULL
    if (is.numeric(y))
      t.ret <- min(y, na.rm = TRUE)
    t.ret
  })

  all.max <- lapply(frame, function(y) {
    t.ret <- NULL
    if (is.numeric(y))
      t.ret <- max(y, na.rm = TRUE)
    t.ret
  })

  list(iord = iord,
       iuno = iuno,
       icon = icon,
       inumord = inumord,
       all.lev = all.lev,
       all.ulev = all.ulev,
       all.dlev = all.dlev,
       all.nlev = all.nlev,
       all.min = all.min,
       all.max = all.max)
}

npKernelBoundsResolve <- function(dati,
                                  varnames,
                                  kerbound = c("none", "range", "fixed"),
                                  kerlb = NULL,
                                  kerub = NULL,
                                  argprefix = "cker") {
  kerbound <- match.arg(kerbound)
  icon.idx <- which(dati$icon)
  ncon <- length(icon.idx)
  if (is.null(varnames) || length(varnames) != length(dati$icon))
    varnames <- paste0("V", seq_along(dati$icon))

  full.lb <- rep(-Inf, length(dati$icon))
  full.ub <- rep(Inf, length(dati$icon))

  if (ncon == 0L) {
    return(list(bound = "none", lb = full.lb, ub = full.ub))
  }

  mins <- unlist(dati$all.min[icon.idx], use.names = FALSE)
  maxs <- unlist(dati$all.max[icon.idx], use.names = FALSE)
  cnames <- varnames[icon.idx]

  recycleBounds <- function(x, nm) {
    if (is.null(x))
      stop(sprintf("'%s' requires '%s' and '%s'.",
                   paste0(argprefix, "bound"),
                   paste0(argprefix, "lb"),
                   paste0(argprefix, "ub")))
    x <- as.numeric(x)
    if (!(length(x) %in% c(1L, ncon))) {
      stop(sprintf("length(%s) must be 1 or equal to the number of continuous variables (%d).",
                   nm, ncon))
    }
    if (length(x) == 1L)
      rep(x, ncon)
    else
      x
  }

  if (kerbound == "none") {
    lb <- rep(-Inf, ncon)
    ub <- rep(Inf, ncon)
  } else if (kerbound == "range") {
    lb <- mins
    ub <- maxs
  } else {
    lb <- recycleBounds(kerlb, paste0(argprefix, "lb"))
    ub <- recycleBounds(kerub, paste0(argprefix, "ub"))
  }

  if (any(!is.finite(lb) & !is.infinite(lb)) || any(!is.finite(ub) & !is.infinite(ub))) {
    stop(sprintf("'%s' and '%s' must be finite values or +/-Inf.",
                 paste0(argprefix, "lb"),
                 paste0(argprefix, "ub")))
  }

  if (any(lb >= ub)) {
    stop(sprintf("Invalid bounds for '%s': require lower < upper for every continuous variable.",
                 paste0(argprefix, "bound")))
  }

  bad.cover <- (lb > mins) | (ub < maxs)
  if (any(bad.cover)) {
    bad.vars <- paste(cnames[bad.cover], collapse = ", ")
    stop(sprintf("Invalid bounds for '%s': require lower <= min(training) and upper >= max(training) for each continuous variable. Violations: %s",
                 paste0(argprefix, "bound"), bad.vars))
  }

  if (kerbound == "fixed" && all(is.infinite(lb) & (lb < 0)) && all(is.infinite(ub) & (ub > 0))) {
    kerbound <- "none"
  }

  full.lb[icon.idx] <- lb
  full.ub[icon.idx] <- ub

  list(bound = kerbound, lb = full.lb, ub = full.ub)
}

npKernelBoundsCheckEval <- function(evaldat,
                                    icon,
                                    kerlb,
                                    kerub,
                                    argprefix = "cker") {
  if (is.null(evaldat) || is.null(icon) || !any(icon))
    return(invisible(TRUE))

  if (is.null(kerlb) || is.null(kerub))
    return(invisible(TRUE))

  icon.idx <- which(icon)
  ncon <- length(icon.idx)
  if (ncon == 0L)
    return(invisible(TRUE))

  lb <- as.numeric(kerlb[icon.idx])
  ub <- as.numeric(kerub[icon.idx])

  if (all(is.infinite(lb) & lb < 0) && all(is.infinite(ub) & ub > 0))
    return(invisible(TRUE))

  enames <- names(evaldat)
  if (is.null(enames) || length(enames) != length(icon))
    enames <- paste0("V", seq_along(icon))
  cnames <- enames[icon.idx]

  econ <- as.matrix(evaldat[, icon.idx, drop = FALSE])
  if (!is.double(econ))
    storage.mode(econ) <- "double"

  bad.low <- vapply(seq_len(ncon), function(i) {
    is.finite(lb[i]) && any(econ[, i] < lb[i], na.rm = TRUE)
  }, logical(1L))
  bad.high <- vapply(seq_len(ncon), function(i) {
    is.finite(ub[i]) && any(econ[, i] > ub[i], na.rm = TRUE)
  }, logical(1L))

  if (any(bad.low) || any(bad.high)) {
    low.msg <- if (any(bad.low)) {
      paste(sprintf("%s < %s", cnames[bad.low], format(lb[bad.low], digits = 12)), collapse = ", ")
    } else {
      NULL
    }
    high.msg <- if (any(bad.high)) {
      paste(sprintf("%s > %s", cnames[bad.high], format(ub[bad.high], digits = 12)), collapse = ", ")
    } else {
      NULL
    }
    msg <- paste(c(low.msg, high.msg), collapse = "; ")
    stop(sprintf("Evaluation data violate '%s' bounds: %s",
                 paste0(argprefix, "bound"), msg))
  }

  invisible(TRUE)
}

npKernelBoundsMarshal <- function(kerlb, kerub) {
  lb <- as.double(kerlb)
  ub <- as.double(kerub)
  big <- .Machine$double.xmax
  lb[is.infinite(lb) & lb < 0] <- -big
  ub[is.infinite(ub) & ub > 0] <- big
  list(lb = lb, ub = ub)
}

validateBandwidth <- function(bws){
  vari <- names(bws$bandwidth)
  bchecker <- function(j){
    v <- vari[j]
    dati <- bws$dati[[v]]
    bwv <- bws$bandwidth[[j]]
    stopifnot(length(bwv) == length(dati$iord))

    cd <- function(a,b){
      (a-b)/(a+b+.Machine$double.eps) > 5.0*.Machine$double.eps
    }
    
    vb <- sapply(1:length(bwv), function(i){
      falg <- (bwv[i] < 0)

      if (dati$icon[i] && (falg || (!is.finite(bwv[i])))){
        stop(paste("Invalid bandwidth supplied for continuous",
                   "variable", bws$varnames[[v]][i], ":",bwv[i]))
      }
      
      if (dati$iord[i] &&
          (falg || cd(bwv[i],oMaxL(dati$all.nlev[[i]],
                         kertype = bws$klist[[v]]$okertype)))){
        stop(paste("Invalid bandwidth supplied for ordered",
                   "variable", bws$varnames[[v]][i], ":",bwv[i]))
      }
      
      if (dati$iuno[i] &&
          (falg || cd(bwv[i],uMaxL(dati$all.nlev[[i]],
                         kertype = bws$klist[[v]]$ukertype)))){
        stop(paste("Invalid bandwidth supplied for unordered",
                   "variable", bws$varnames[[v]][i], ":",bwv[i]))
      }
      return(TRUE)
    })
    
    return(vb)
  }
  vbl <- lapply(1:length(vari), bchecker)
  invisible(vbl)
}

validateBandwidthTF <- function(bws){
  vari <- names(bws$bandwidth)
  bchecker <- function(j){
    v <- vari[j]
    dati <- bws$dati[[v]]
    bwv <- bws$bandwidth[[j]]

    if(length(bwv) != length(dati$iord))
      return(FALSE)

    cd <- function(a,b){
      (a-b)/(a+b+.Machine$double.eps) > 5.0*.Machine$double.eps
    }
    
    vb <- sapply(1:length(bwv), function(i){
      falg <- (bwv[i] < 0)

      if (dati$icon[i]) {
        if(bws$type == "fixed") {
          if(falg || (!is.finite(bwv[i]))){
            return(FALSE)
          }
        } else if((bwv[i] < 1) || (!is.finite(bwv[i]))) {
          return(FALSE)
        }
      }
      
      if (dati$iord[i] &&
          (falg || cd(bwv[i],oMaxL(dati$all.nlev[[i]],
                         kertype = bws$klist[[v]]$okertype)))){
        return(FALSE)
      }
      
      if (dati$iuno[i] &&
          (falg || cd(bwv[i],uMaxL(dati$all.nlev[[i]],
                         kertype = bws$klist[[v]]$ukertype)))){
        return(FALSE)
      }
      return(TRUE)
    })
    
    return(all(vb))
  }
  vbl <- all(unlist(lapply(1:length(vari), bchecker)))
  return(vbl)
}


explodeFormula <- function(formula, data=NULL){
  if(any(grepl("\\.",deparse(formula)))) {
      if(is.null(data)) stop("'.' in formula and no 'data' argument")
      formula <- terms(formula, data=data)
  }
  res <- strsplit(strsplit(paste(deparse(formula), collapse=""),
                           " *[~] *")[[1]], " *[+] *")
  stopifnot(all(sapply(res,length) > 0))
  names(res) <- c("response","terms")
  res
}


explodePipe <- function(formula, env = parent.frame()){
  if (!inherits(formula, "formula"))
    formula <- eval(formula, env)
  tf <- as.character(formula)  
  tf <- tf[length(tf)]
  lhs <- if (length(as.character(formula)) == 3) {
    strsplit(as.character(formula)[2], " *[+] *")
  } else {
    list()
  }
  rhs <- strsplit(strsplit(tf, " *[|] *")[[1]], " *[+] *")
  c(lhs, rhs)
}

"%~%" <- function(a,b) {
  identical(class(a), class(b)) && (length(a) == length(b)) &&
  all(unlist(lapply(a,coarseclass)) == unlist(lapply(b,coarseclass)))
}

coarseclass <- function(a) {
  if (inherits(a, "integer")) return("numeric")
  return(class(a)[1])
}

toFrame <- function(frame) {
  if(!is.data.frame(frame)){
    t.names <- NULL

    if(!(is.vector(frame) || is.factor(frame) || is.matrix(frame)))
      stop(deparse(substitute(frame))," must be a data frame, matrix, vector, or factor")

    if(!is.matrix(frame))
      t.names <- deparse(eval(substitute(substitute(frame)), envir = parent.frame()))
    
    frame <- data.frame(frame, check.names=FALSE)
    
    if(!is.null(t.names))
      names(frame) <- t.names
  }
  return(frame)
}


cast <- function(a, b, same.levels = TRUE){
  if(is.ordered(b)){
    if(same.levels)
      ordered(a, levels = levels(b))
    else
      ordered(a)
  }   
  else if(is.factor(b)){
    if(same.levels)
      factor(a, levels = levels(b))
    else
      factor(a)
  }
  else if (coarseclass(b) == "numeric")
    as.double(a)
  else if (is.data.frame(b)) {
    if (dim(a)[2] == dim(b)[2]){
      r = data.frame(a)
      for (i in 1: length(b))
        r[,i] = cast(a[,i],b[,i], same.levels = same.levels)
      r
    } else { stop("a could not be cast as b") }
  }
}

subcol <- function(x, v, i){
  x[,i] = cast(v,x[,i])
  x
}

mcvConstruct <- function(dati){
  nuno <- sum(dati$iuno)
  nord <- sum(dati$iord)

  num.row <- max(sapply(dati$all.lev,length))

  pad.num <- numNotIn(unlist(dati$all.dlev))

  mcv <- matrix(data = pad.num, nrow = num.row, ncol = (nuno+nord))
  attr(mcv, "num.row") <- num.row
  attr(mcv, "pad.num") <- pad.num

  cnt <- 0
  if (nuno > 0)
    for (i in which(dati$iuno)) 
      mcv[1:length(dati$all.lev[[i]]), (cnt <- cnt+1)] <- dati$all.dlev[[i]]

  cnt <- 0
  if (nord > 0)
    for (i in which(dati$iord))
      mcv[1:length(dati$all.lev[[i]]), (cnt <- cnt+1)+nuno] <- dati$all.dlev[[i]]

  mcv
}

## when admitting new categories, adjustLevels will attempt to catch possible mistakes:
## if an unordered variable contains more than one new category, warn
## if an ordered, but scaleless variable contains a new category, error
## if an ordered, scale-possessing variable contains a new category lacking scale, error

adjustLevels <- function(data, dati, allowNewCells = FALSE){
  for (i in which(dati$iord | dati$iuno)){
    if (allowNewCells){
      newCats <- setdiff(levels(data[,i]),dati$all.lev[[i]])
      if (length(newCats) >= 1){
        if (dati$iuno[i]){
          if (length(newCats) > 1)
            warning(paste("more than one 'new' category is redundant when estimating on unordered data.\n",
                          "training data categories: ", paste(dati$all.lev[[i]], collapse=" "),"\n",
                          "redundant estimation data categories: ", paste(newCats, collapse=" "), "\n", sep=""))
          data[,i] <- factor(data[,i], levels = c(dati$all.lev[[i]], newCats))
        } else {
          if (dati$inumord[i]){
            if (!isNum(newCats))
              stop(paste("estimation data contains a new qualitative category, but training data is\n",
                         "ordered, and numeric.\n",
                         "training data categories: ", paste(dati$all.lev[[i]], collapse=" "),"\n",
                         "conflicting estimation data categories: ", paste(newCats, collapse=" "), "\n", sep=""))
          } else {
            stop(paste("estimation beyond the support of training data of an ordered,\n",
                       "categorical, qualitative variable is not supported.\n"))
          }

          data[,i] <- ordered(data[,i], levels = sort(as.numeric(c(dati$all.lev[[i]], newCats))))
        }
      } else {
        data[,i] <- factor(data[,i], levels = dati$all.lev[[i]])
      }
    } else {
      if (!all(is.element(levels(data[,i]), dati$all.lev[[i]])))
        stop("data contains unknown factors (wrong dataset provided?)")
      data[,i] <- factor(data[,i], levels = dati$all.lev[[i]])
    }
  }

  data
}

toMatrix <- function(data) {
  tq <- sapply(data, function(y) {
    if(is.factor(y))
      y <- dlev(y)[as.integer(y)]
    y})
  dim(tq) <- dim(data) ## cover the corner case of single element d.f.
  tq
}

## this doesn't just strictly check for the response, but does tell you
## that evaluating with response fails ... in principle the evaluation
## could fail without the response too, but then the calling routine is about
## to die a noisy death anyhow ...
succeedWithResponse <- function(tt, frame){
  !inherits(try(eval(expr = attr(tt, "variables"),
                     envir = frame, enclos = NULL), silent = TRUE), "try-error")
}

## determine whether a bandwidth
## matches a data set
bwMatch <- function(data, dati){
  if (length(dati$icon) != ncol(data))
    stop("bandwidth vector is of improper length")

  test.dati <- untangle(data)

  if (any(xor(dati$icon,test.dati$icon)) ||
      any(xor(dati$iord,test.dati$iord)) ||
      any(xor(dati$iuno,test.dati$iuno)))
    stop(paste("supplied bandwidths do not match","data", "in type"))
}

uMaxL <- function(c, kertype = c("aitchisonaitken","liracine")){
  switch(kertype,
         aitchisonaitken = (c-1.0)/c,
         liracine = 1.0)
}

oMaxL <- function(c, kertype = c("wangvanryzin", "liracine", "nliracine", "racineliyan")){
  switch(kertype,
         wangvanryzin = 1.0,
         liracine = 1.0,
         nliracine = 1.0,
         "racineliyan" = 1.0)
}

## tested with: rbandwidth
## right now all bandwidth objects have some crusty
## vestiges of their evolution, ie. non-list metadata
## such as xnames or ynames. The new metadata system is
## for the most part list based and facilitates generic
## operations

updateBwNameMetadata <- function(nameList, bws){
  ## names of 'old' metadata in bw object
  onames <- names(nameList)
  lapply(1:length(nameList), function(i) {
    bws[[onames[i]]] <<- nameList[[i]]
    bws$varnames[[substr(onames[i],1,1)]] <<- nameList[[i]]
  })
  return(bws)
}

## some string utility functions

pad <- function(s){
  ifelse(nchar(s) > 0, paste("",s,""), " ")
}

rpad <- function(s){
  ifelse(nchar(s) > 0, paste(s,""), "")
}

blank <- function(len){
  sapply(len, function(nb){
    paste(rep(' ', times = nb), collapse='')
  })
}

formatv <- function(v){
  sapply(1:length(v), function(j){ format(v[j]) })
}

## strings used in various report generating functions

genOmitStr <- function(x){
  t.str <- ''
  if(!is.null(x$rows.omit) & !identical(x$rows.omit, NA))
    t.str <- paste("\nNo. Complete Observations: ", x$nobs,
                   "\nNo. Incomplete (NA) Observations: ", x$nobs.omit,
                   "\nObservations omitted or excluded: ", paste(x$rows.omit, collapse=" "),
                   "\n")
  return(t.str)
}

## Estimation-related rgf's
genGofStr <- function(x){
  paste(ifelse(is.na(x$MSE),"",paste("\nResidual standard error:",
                                     format(sqrt(x$MSE)))),
        ifelse(is.na(x$R2),"",paste("\nR-squared:",
                                    format(x$R2))), sep="")
}

genTimingStr <- function(x){
  if (is.null(x$total.time) || is.na(x$total.time))
    return("")

  if (!is.null(x$optim.time) && !is.na(x$optim.time) &&
      !is.null(x$fit.time) && !is.na(x$fit.time))
    return(paste("\nEstimation Time: ", format(x$total.time), " seconds (optim ",
                 format(x$optim.time), "s, fit ", format(x$fit.time), "s)", sep = ""))

  paste("\nEstimation Time: ",format(x$total.time)," seconds",sep = "")
}
  
pCatGofStr <- function(x){
  if(!identical(x$confusion.matrix, NA)){
    cat("\nConfusion Matrix\n")
    print(x$confusion.matrix)
  }

  if (!identical(x$CCR.overall,NA))
    cat("\nOverall Correct Classification Ratio: ", format(x$CCR.overall))

  if (!identical(x$CCR.byoutcome,NA)){
    cat("\nCorrect Classification Ratio By Outcome:\n")
    print(x$CCR.byoutcome)
  }

  if (!identical(x$fit.mcfadden,NA))
    cat("\nMcFadden-Puig-Kerschner performance measure: ", format(x$fit.mcfadden))

}

genDenEstStr <- function(x){
  paste("\nBandwidth Type: ",x$ptype,
        ifelse(is.null(x$log_likelihood) || identical(x$log_likelihood, NA),"",
               paste("\nLog Likelihood:",
                     format(x$log_likelihood))),
        sep="")
}

genRegEstStr <- function(x){
  regtype <- if (!is.null(x$regtype)) x$regtype else if (!is.null(x$bws)) x$bws$regtype else NULL
  basis <- if (!is.null(x$basis)) x$basis else if (!is.null(x$bws)) x$bws$basis else NULL
  bern <- if (!is.null(x$bernstein.basis)) x$bernstein.basis else if (!is.null(x$bws)) x$bws$bernstein.basis else NULL
  est.label <- if (identical(regtype, "lp")) npFormatRegressionType(x) else x$pregtype
  basis.family <- if (identical(regtype, "lp")) npLpBasisFamilyLabel(basis) else NULL
  basis.rep <- if (identical(regtype, "lp")) npLpBasisRepresentationLabel(bern) else NULL
  paste(ifelse(is.null(est.label),"",
               paste("\nKernel Regression Estimator:",est.label)),
        ifelse(is.null(basis.family), "",
               paste("\nLP Basis Family:", basis.family)),
        ifelse(is.null(basis.rep), "",
               paste("\nLP Basis Representation:", basis.rep)),
        ifelse(is.null(x$ptype), "",
               paste("\nBandwidth Type:",x$ptype)),
        ifelse(is.null(x$tau), "", paste("\nTau:", x$tau)),
        sep = "")
}

npLpBasisFamilyLabel <- function(basis){
  b <- tolower(ifelse(is.null(basis) || !length(basis), "glp", basis))
  switch(b,
         glp = "Generalized",
         additive = "Additive",
         tensor = "Tensor",
         tools::toTitleCase(b))
}

npLpBasisRepresentationLabel <- function(bernstein){
  if (isTRUE(bernstein)) "Bernstein" else "Raw"
}

npLpBasisNcol <- function(basis = "glp", degree){
  if (is.null(degree) || !length(degree))
    return(NA_real_)
  d <- dim_basis(basis = basis,
                 kernel = TRUE,
                 degree = as.integer(degree),
                 segments = rep.int(1L, length(degree)))
  if (identical(tolower(basis), "tensor")) d else d + 1.0
}

npFormatRegressionType <- function(x){
  regtype <- if (!is.null(x$regtype)) {
    x$regtype
  } else if (!is.null(x$bws) && !is.null(x$bws$regtype)) {
    x$bws$regtype
  } else {
    NULL
  }

  pregtype <- if (!is.null(x$pregtype)) {
    x$pregtype
  } else if (!is.null(x$bws) && !is.null(x$bws$pregtype)) {
    x$bws$pregtype
  } else {
    NULL
  }

  if (!identical(regtype, "lp"))
    return(pregtype)

  degree <- if (!is.null(x$degree)) {
    x$degree
  } else if (!is.null(x$bws) && !is.null(x$bws$degree)) {
    x$bws$degree
  } else {
    NULL
  }

  if (is.null(degree) || length(degree) == 0)
    return("Local-Polynomial")

  basis <- if (!is.null(x$basis)) {
    x$basis
  } else if (!is.null(x$bws) && !is.null(x$bws$basis)) {
    x$bws$basis
  } else {
    "glp"
  }
  basis.family <- npLpBasisFamilyLabel(basis)

  sprintf("Local-Polynomial (%s basis; degree = %s)",
          basis.family, paste(degree, collapse = ","))
}


## bandwidth-related report generating functions
genBwSelStr <- function(x){
  fval.str <- ""
  if (!identical(x$fval, NA)) {
    fval.str <- if (is.null(x$ifval) || identical(x$ifval, NA)) {
      paste("\nObjective Function Value: ", format(x$fval), sep = "")
    } else {
      paste("\nObjective Function Value: ", format(x$fval),
            " (achieved on multistart ", x$ifval, ")", sep = "")
    }
  }

  nfe.str <- ""
  if(!(is.null(x$num.feval) || identical(x$num.feval, NA))){
    nfe.str <- paste("\nNumber of Function Evaluations: ", format(x$num.feval), sep="")
    if(!(is.null(x$num.feval.fast) || identical(x$num.feval.fast, NA)) &&
       !(is.null(x$num.feval.fallback) || identical(x$num.feval.fallback, NA))){
      nfe.str <- paste(nfe.str,
                       " (fast = ", format(x$num.feval.fast),
                       ", fallback = ", format(x$num.feval.fallback), ")",
                       sep="")
    } else if(!(is.null(x$num.feval.fast) || identical(x$num.feval.fast, NA))){
      nfe.str <- paste(nfe.str, " (fast = ", format(x$num.feval.fast), ")", sep="")
    } else if(!(is.null(x$num.feval.fallback) || identical(x$num.feval.fallback, NA))){
      nfe.str <- paste(nfe.str, " (fallback = ", format(x$num.feval.fallback), ")", sep="")
    }
  }

  pregtype <- npFormatRegressionType(x)

  paste(ifelse(is.null(pregtype),"",paste("\nRegression Type:", pregtype)),
        ifelse(is.null(x$pmethod),"",paste("\nBandwidth Selection Method:",
                                           x$pmethod)),
        if (!identical(x$formula,NULL)) paste("\nFormula:",
                                              paste(deparse(x$formula), collapse="\n")),
        ifelse(is.null(x$ptype), "",
               paste("\nBandwidth Type: ",x$ptype, sep="")),
        fval.str,
        nfe.str,
        sep="")
}

genBwScaleStrs <- function(x){
  ## approach is to take metadata and flatten it so it then can be
  ## processed into a single string

  vari <- names(x$sumNum)
  
  t.icon <- lapply(vari, function(v){
    x$dati[[v]]$icon })

  sumText <- lapply(1:length(vari), function(i) {
    ifelse(t.icon[[i]],
           ifelse(x$type == "fixed",
                  "Scale Factor:",""), "Lambda Max:")
    
  })

  maxNameLen <- max(nchar(unlist(sumText)))
  print.sumText <- lapply(sumText, '!=', "")

  sumText <- lapply(1:length(sumText), function(i){
    paste(blank(maxNameLen - nchar(sumText[[i]])), sumText[[i]], sep="")
  })

  t.nchar <- lapply(x$varnames[vari], nchar)
                    
  maxNameLen <- max(unlist(t.nchar))

  vatText <- lapply(1:length(t.nchar), function(j){
    paste("\n", rpad(x$vartitleabb[[vari[j]]]), "Var. Name: ",
          x$varnames[[vari[j]]],
          sapply(t.nchar[[j]], function(nc){
            paste(rep(' ', maxNameLen - nc), collapse='')
          }), sep='')
  })

  return(sapply(1:length(t.nchar), function(j){
    paste(vatText[[j]], " Bandwidth: ", npFormat(x$bandwidth[[j]]), " ",
          ifelse(print.sumText[[j]],
                 paste(sumText[[j]], " ", npFormat(x$sumNum[[j]]), sep=""), ""),
          sep="", collapse="")
  }))
}

npFormat <- function(x){
  format(sapply(x,format))
}

genBwKerStrs <- function(x){
  vari <- names(x$klist)

  ncon <- sapply(vari, function(v){
    sum(x$dati[[v]]$icon)
  })

  nuno <- sapply(vari, function(v){
    sum(x$dati[[v]]$iuno)
  })

  nord <- sapply(vari, function(v){
    sum(x$dati[[v]]$iord)
  })

  cktype <- sapply(vari, function(v){
    x$klist[[v]]$ckertype
  })

  uktype <- sapply(vari, function(v){
    x$klist[[v]]$ukertype
  })

  oktype <- sapply(vari, function(v){
    x$klist[[v]]$okertype
  })

  tt <- ''

  if(any(ncon > 0)){
    tt <- paste("\n",
                ifelse(length(unique(cktype)) == 1,
                       paste("\nContinuous Kernel Type:",
                             x$klist[[vari[1]]]$pckertype),
                       paste(sapply(1:length(vari), function(v){
                         ifelse(ncon[v] > 0,
                                paste("\nContinuous Kernel Type (",
                                      x$vartitleabb[[vari[v]]],
                                      " Var.): ", x$klist[[vari[v]]]$pckertype, sep=""),"")
                       }), collapse = "")),
                sep = "")
    tt <-
      paste(tt, paste(sapply(1:length(vari), function(i){
        ifelse(ncon[i] > 0,
               paste("\nNo. Continuous", pad(x$vartitle[[vari[i]]]), "Vars.: ",
                     ncon[i], sep=""), "")
      }), collapse = ""), sep="")
  }
                
    
  if(any(nuno > 0)) {
    tt <- paste(tt, "\n",
                ifelse(length(unique(uktype)) == 1,
                       paste("\nUnordered Categorical Kernel Type:",
                             x$klist[[vari[1]]]$pukertype),
                       paste(sapply(1:length(vari), function(i){
                         ifelse(nuno[i] > 0,
                                paste("\nUnordered Categorical Kernel Type (",
                                      x$vartitleabb[[vari[i]]],
                                      " Var.): ", x$klist[[vari[i]]]$pukertype, sep=""),"")
                       }), collapse = "")),
                sep = "")
    tt <-
      paste(tt, paste(sapply(1:length(vari), function(i){
        ifelse(nuno[i] > 0,
               paste("\nNo. Unordered Categorical", pad(x$vartitle[[vari[i]]]), "Vars.: ",
                     nuno[i], sep=""), "")
      }), collapse = ""), sep="")

  }

  if(any(nord > 0)) {
    tt <- paste(tt, "\n",
                ifelse(length(unique(oktype)) == 1,
                paste("\nOrdered Categorical Kernel Type:",
                      x$klist[[vari[1]]]$pokertype),
                paste(sapply(1:length(vari), function(i){
                  ifelse(nord[i] > 0,
                         paste("\nOrdered Categorical Kernel Type (",
                               x$vartitleabb[[vari[i]]],
                               " Var.): ", x$klist[[vari[i]]]$pokertype, sep=""),"")
                }), collapse = "")),
         sep = "")
    tt <-
      paste(tt, paste(sapply(1:length(vari), function(i){
        ifelse(nord[i] > 0,
               paste("\nNo. Ordered Categorical", pad(x$vartitle[[vari[i]]]), "Vars.: ",
                     nord[i], sep=""), "")
      }), collapse = ""), sep="")

  }

  return(tt)
}

genBwKerStrsXY <- function(x){
  t.str <- ''
  cnt <- 0
  
  if (x$xncon + x$yncon > 0){
    t.str[cnt <- cnt + 1] <-
      paste(ifelse(x$pcxkertype == x$pcykertype,
                   paste("\n\nContinuous Kernel Type:",x$pcxkertype),
                   paste("\n", ifelse(x$xncon > 0,
                                      paste("\nContinuous Kernel Type (Exp. Var.):",
                                            x$pcxkertype), ""),
                         ifelse(x$yncon > 0,
                                paste("\nContinuous Kernel Type (Dep. Var.):",
                                      x$pcykertype), ""))),
            ifelse(x$yncon > 0,paste("\nNo. Continuous Dependent Vars.:",x$yncon),""),
            ifelse(x$xncon > 0,paste("\nNo. Continuous Explanatory Vars.:",x$xncon),""))
  }

  if (x$xnuno + x$ynuno > 0){
    t.str[cnt <- cnt + 1] <-
      paste(ifelse(x$puxkertype == x$puykertype,
                   paste("\n\nUnordered Categorical Kernel Type:",x$puxkertype),
                   paste("\n", ifelse(x$xnuno > 0,
                                      paste("\nUnordered Categorical Kernel Type (Exp. Var.):",
                                            x$puxkertype), ""),
                         ifelse(x$ynuno > 0,
                                paste("\nUnordered Categorical Kernel Type (Dep. Var.):",
                                      x$puykertype), ""))),
            ifelse(x$ynuno > 0,paste("\nNo. Unordered Categorical Dependent Vars.:",x$ynuno),""),
            ifelse(x$xnuno > 0,paste("\nNo. Unordered Categorical Explanatory Vars.:",x$xnuno),""))
  }

  if (x$xnord + x$ynord > 0){
    t.str[cnt <- cnt + 1] <-
      paste(ifelse(x$poxkertype == x$poykertype,
                   paste("\n\nOrdered Categorical Kernel Type:",x$poxkertype),
                   paste("\n", ifelse(x$xnord > 0,
                                      paste("\nOrdered Categorical Kernel Type (Exp. Var.):",
                                            x$poxkertype), ""),
                         ifelse(x$ynord > 0,
                                paste("\nOrdered Categorical Kernel Type (Dep. Var.):",
                                      x$poykertype), ""))),
            ifelse(x$ynord > 0,paste("\nNo. Ordered Categorical Dependent Vars.:",x$ynord),""),
            ifelse(x$xnord > 0,paste("\nNo. Ordered Categorical Explanatory Vars.:",x$xnord),""))
  }
  return(t.str)
}

genBwGOFStrs <- function(x) {
  ###paste("Residual standard error:", sqrt
}
  
## statistical functions
RSQfunc <- function(y,y.pred) {
  y.mean <- mean(y)
  (sum((y-y.mean)*(y.pred-y.mean))^2)/(sum((y-y.mean)^2)*sum((y.pred-y.mean)^2))
}

MSEfunc <- function(y,y.fit) {
  mean((y-y.fit)^2)
}

MAEfunc <- function(y,y.fit) {
  mean(abs(y-y.fit))
}

MAPEfunc <- function(y,y.fit) {
  jj = which(y != 0)
  
  mean(c(abs((y[jj]-y.fit[jj])/y[jj]), as.numeric(replicate(length(y)-length(jj),2))))
}

CORRfunc <- function(y,y.fit) {
  abs(corr(cbind(y,y.fit)))
}

SIGNfunc <- function(y,y.fit) {
  sum(sign(y) == sign(y.fit))/length(y)
}

EssDee <- function(y){
  if(any(dim(as.matrix(y)) == 0))
      return(0)
  sd.vec <- apply(as.matrix(y),2,sd)
  IQR.vec <- apply(as.matrix(y),2,IQR)/QFAC
  mad.vec <- apply(as.matrix(y),2,mad)
  a <- apply(cbind(sd.vec,IQR.vec,mad.vec),1, function(x) max(x))
  if(any(a<=0)) warning(paste("variable ",which(a<=0)," appears to be constant",sep=""))
  a <- apply(cbind(sd.vec,IQR.vec,mad.vec),1, function(x) min(x[x>0]))  
  return(a)
}


##EssDee <- function(y){
##
##  sd.vec <- apply(as.matrix(y),2,sd)
##  IQR.vec <- apply(as.matrix(y),2,IQR)/QFAC
##  return(ifelse(sd.vec<IQR.vec|IQR.vec==0,sd.vec,IQR.vec))
##  
##}

## consolidating various bits of code related to converting internal settings
## to printable strings

bwmToPrint <- function(s){
  switch(s,
         manual = "Manual",
         cv.aic = "Expected Kullback-Leibler Cross-Validation",
         cv.ml = "Maximum Likelihood Cross-Validation",
         cv.cdf = "Least Squares Cross-Validation",
         cv.ls = "Least Squares Cross-Validation",
         cv.ls.np = "Least Squares Cross-Validation (block algorithm)",
         "normal-reference" = "Normal Reference")
}

bwtToPrint <- function(s){
  switch(s,
         fixed = "Fixed",
         generalized_nn = "Generalized Nearest Neighbour",
         adaptive_nn = "Adaptive Nearest Neighbour" )
}

cktToPrint <- function(s, order = "", kerbound = "none"){
  pck <- switch(s,
                gaussian = paste(order,"Gaussian"),
                epanechnikov =  paste(order,"Epanechnikov"),
                uniform = "Uniform",
                "truncated gaussian" = "Truncated Gaussian")
  if (!is.null(kerbound) && !identical(kerbound, "none"))
    pck <- paste0(pck, " (bounded/", kerbound, ")")
  pck
}

uktToPrint <- function(s){
  switch(s,
         aitchisonaitken = "Aitchison and Aitken",
         liracine = "Li and Racine (normalized)")
}

oktToPrint <- function(s, normalized = FALSE) {
  if(normalized){
    pok <- 
      switch(s,
             wangvanryzin = "Wang and Van Ryzin", 
             liracine = "Li and Racine (normalized)",
             nliracine = "Li and Racine (normalized)",
             "racineliyan" = "Racine, Li, and Yan")
  } else {
    pok <- 
      switch(s,
             wangvanryzin = "Wang and Van Ryzin", 
             liracine = "Li and Racine",
             nliracine = "Li and Racine (normalized)",
             "racineliyan" = "Racine, Li, and Yan")
  }
  return(pok)
}

### holding place for some generic methods

se <- function(x){
  UseMethod("se",x)
}

gradients <- function(x, ...){
  UseMethod("gradients",x)
}

## From crs to avoid crs:::W.glp

mypoly <- function(x,
                   ex=NULL,
                   degree,
                   gradient.compute = FALSE,
                   r=0,
                   Bernstein = TRUE) {

  if(missing(x)) stop(" Error: x required")
  if(missing(degree)) stop(" Error: degree required")
  if(degree < 1) stop(" Error: degree must be a positive integer")
  if(!is.logical(Bernstein)) stop(" Error: Bernstein must be logical")

  if(!Bernstein) {

    ## Raw polynomials and their derivatives

    if(!is.null(ex)) x <- ex

    if(gradient.compute) {
      Z <- NULL
      for(i in 1:degree) {
        if((i-r) >= 0) {
          tmp <- (factorial(i)/factorial(i-r))*x^max(0,i-r)
        } else {
          tmp <- rep(0,length(x))
        }
        Z <- cbind(Z,tmp)
      }
    } else {
      Z <- outer(x,1L:degree,"^")
    }

  } else {
    ## Bernstein polynomials and their derivatives (i.e. Bezier curves
    ## i.e. B-splines with no interior knots)
    if(is.null(ex)) {
      if(gradient.compute) {
        Z <- gsl.bs(x,degree=degree,deriv=r)
      } else {
        Z <- gsl.bs(x,degree=degree)
      }
    } else {
      if(gradient.compute) {
        Z <- predict(gsl.bs(x,degree=degree,deriv=r),newx=ex)
      } else {
        Z <- predict(gsl.bs(x,degree=degree),newx=ex)
      }
    }
  }

  return(as.matrix(Z))

}

## W.lp accepts a vector of degrees and provides local-polynomial bases
## with selectable term structure (glp/additive/tensor).

npBuildLpTerms <- function(degree, basis = c("glp", "additive", "tensor")) {
  basis <- match.arg(basis)
  k <- length(degree)
  if (k == 0L)
    return(matrix(integer(0), nrow = 1L, ncol = 0L))

  degree <- as.integer(degree)
  degree.list <- lapply(degree, function(d) 0:d)
  z <- as.matrix(do.call("expand.grid", degree.list))
  s <- rowSums(z)

  if (identical(basis, "glp")) {
    ind <- (s > 0) & (s <= max(degree))
    z <- z[ind, , drop = FALSE]
    if (!all(degree == max(degree))) {
      for (j in seq_along(degree)) {
        d <- degree[j]
        if ((d < max(degree)) && (d > 0)) {
          s <- rowSums(z)
          dropj <- (s > d) & (z[, j, drop = FALSE] == matrix(d, nrow(z), 1, byrow = TRUE))
          z <- z[!dropj, , drop = FALSE]
        }
      }
    }
  } else if (identical(basis, "additive")) {
    ind <- (s > 0) & (rowSums(z > 0) == 1L)
    z <- z[ind, , drop = FALSE]
  } else if (identical(basis, "tensor")) {
    z <- z[s > 0, , drop = FALSE]
  }

  rbind(matrix(0L, nrow = 1L, ncol = k), z)
}

W.lp <- function(xdat = NULL,
                 exdat = NULL,
                 degree = NULL,
                 gradient.vec = NULL,
                 basis = c("glp", "additive", "tensor"),
                 bernstein.basis = TRUE,
                 Bernstein = bernstein.basis) {

  if(is.null(xdat)) stop(" Error: You must provide data")
  if(is.null(degree) || any(degree < 0)) stop(paste(" Error: degree vector must contain non-negative integers\ndegree is (", degree, ")\n",sep=""))
  basis <- match.arg(basis)

  xdat <- as.data.frame(xdat)

  xdat.col.numeric <- sapply(1:ncol(xdat),function(i){is.numeric(xdat[,i])})
  k <- ncol(as.data.frame(xdat[,xdat.col.numeric]))

  xdat.numeric <- NULL
  if(k > 0) {
    xdat.numeric <- as.data.frame(xdat[,xdat.col.numeric])
    if(!is.null(exdat)) {
      exdat.numeric <- as.data.frame(exdat[,xdat.col.numeric])
    } else {
      exdat.numeric <- NULL
    }
  }

  if(!is.null(gradient.vec) && (length(gradient.vec) != k)) stop(paste(" Error: gradient vector and number of numeric predictors must be conformable\n",sep=""))
  if(!is.null(gradient.vec) && any(gradient.vec < 0)) stop(paste(" Error: gradient vector must contain non-negative integers\n",sep=""))

  if(!is.null(gradient.vec)) {
    gradient.compute <- TRUE
  } else {
    gradient.compute <- FALSE
    gradient.vec <- rep(NA,k)
  }

  if(length(degree) != k) stop(" Error: degree vector and number of numeric predictors incompatible")

  if(all(degree == 0) || (k == 0)) {

    ## Local constant OR no continuous variables

    if(is.null(exdat)) {
      return(matrix(1,nrow=nrow(xdat.numeric),ncol=1))
    } else {
      return(matrix(1,nrow=nrow(exdat.numeric),ncol=1))
    }

  } else {

    z <- npBuildLpTerms(degree = degree, basis = basis)
    z.noi <- z[-1L, , drop = FALSE]
    if(is.null(exdat)) {
      res <- rep.int(1,nrow(xdat.numeric))
    } else {
      res <- rep.int(1,nrow(exdat.numeric))
    }
    res.deriv <- 1
    if(degree[1] > 0) {
      res <- cbind(1, mypoly(x=xdat.numeric[,1],
                             ex=exdat.numeric[,1],
                             degree=degree[1],
                             gradient.compute=gradient.compute,
                             r=gradient.vec[1],
                             Bernstein=Bernstein))[, 1 + z.noi[, 1]]

      if(gradient.compute && gradient.vec[1] != 0) res.deriv <- cbind(1,matrix(NA,1,degree[1]))[, 1 + z.noi[, 1],drop=FALSE]
      if(gradient.compute && gradient.vec[1] == 0) res.deriv <- cbind(1,matrix(0,1,degree[1]))[, 1 + z.noi[, 1],drop=FALSE]
    }
    if(k > 1) for (i in 2:k) if(degree[i] > 0) {
      res <- res * cbind(1, mypoly(x=xdat.numeric[,i],
                                   ex=exdat.numeric[,i],
                                   degree=degree[i],
                                   gradient.compute=gradient.compute,
                                   r=gradient.vec[i],
                                   Bernstein=Bernstein))[, 1 + z.noi[, i]]
      if(gradient.compute && gradient.vec[i] != 0) res.deriv <- res.deriv * cbind(1,matrix(NA,1,degree[i]))[, 1 + z.noi[, i],drop=FALSE]
      if(gradient.compute && gradient.vec[i] == 0) res.deriv <- res.deriv *cbind(1,matrix(0,1,degree[i]))[, 1 + z.noi[, i],drop=FALSE]
    }

    if(is.null(exdat)) {
      res <- matrix(res,nrow=NROW(xdat))
    } else {
      res <- matrix(res,nrow=NROW(exdat))
    }
    if(gradient.compute) res.deriv <- matrix(res.deriv,nrow=1)
    colnames(res) <- apply(z.noi, 1L, function(x) paste(x, collapse = "."))
    if(gradient.compute) colnames(res.deriv) <- apply(z.noi, 1L, function(x) paste(x, collapse = "."))

    if(gradient.compute) {
      res[,!is.na(as.numeric(res.deriv))] <- 0
      return(cbind(0,res))
    } else {
      return(cbind(1,res))
    }

  }

}

### internal constants used in the c backend

SF_NORMAL = 0
SF_ARB = 1

BW_FIXED = 0
BW_GEN_NN = 1
BW_ADAP_NN = 2

IMULTI_TRUE = 1
IMULTI_FALSE = 0

RE_MIN_TRUE = 0
RE_MIN_FALSE = 1

IO_MIN_TRUE = 1
IO_MIN_FALSE = 0

USE_START_NO = 0
USE_START_YES = 1

NP_DO_DENS = 1
NP_DO_DIST = 0

## initially making an np-wide option via the 'options' mechanism
DO_TREE_NO = 0
DO_TREE_YES = 1

##kernel defs
CKER_GAUSS = 0
CKER_EPAN  = 4
CKER_UNI   = 8
CKER_TGAUSS = 9

UKER_AIT = 0
UKER_LR = 1

OKER_WANG = 0
OKER_LR = 1
OKER_NLR = 2
OKER_RLY = 3

##density 
BWM_CVML = 0
BWM_CVLS = 1
BWM_CVML_NP= 2

##distribution
DBWM_CVLS = 0

##regression
BWM_CVAIC = 0

REGTYPE_LC = 0
REGTYPE_LL = 1
REGTYPE_GLP = 2

##conditional density/distribution
CBWM_CVML = 0
CBWM_CVLS = 1
CBWM_NPLS = 2
CBWM_CCDF = 3 # Added 7/2/2010 jracine

##conditional distribution
CDBWM_CVLS = 0

##integral operators on kernels
OP_NOOP        = -1
OP_NORMAL      = 0
OP_CONVOLUTION = 1
OP_DERIVATIVE  = 2
OP_INTEGRAL    = 3

ALL_OPERATORS = c(OP_NORMAL, OP_CONVOLUTION, OP_DERIVATIVE, OP_INTEGRAL)
names(ALL_OPERATORS) <- c("normal","convolution", "derivative", "integral")

PERMUTATION_OPERATORS <- c(OP_NOOP, OP_NORMAL, OP_DERIVATIVE, OP_INTEGRAL)
names(PERMUTATION_OPERATORS) <- c("none", "normal", "derivative", "integral")

## useful numerical constants of kernel integrals
int.kernels <- c(0.28209479177387814348, 0.47603496111841936711, 0.62396943688265038571, 0.74785078617543927990,
                 0.26832815729997476357, 0.55901699437494742410, 0.84658823667359826246, 1.1329342579014329689,
                 0.5, 2.90113075268188e-01)

QFAC <- qnorm(.25,lower.tail=F)*2
