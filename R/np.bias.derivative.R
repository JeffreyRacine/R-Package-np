# Pure coordinate derivatives of the complete response-linear correction.
# First-order owners are retained; the shared recurrence is used only when
# a newly requested order exceeds one. No training-square operator is formed.
.np_bias_derivative_request <- function(reg, d, order) {
  if (is.null(order)) return(NULL)
  order <- npValidateGlpGradientOrder("lp",order,d)
  if (all(order==1L)) return(NULL)
  if (any(order>1L & order>=reg[["degree.engine"]]))
    stop("higher-order bias correction requires gradient.order <= degree - 1 in each requested higher-order coordinate",call.=FALSE)
  order
}

.np_bias_derivative_order <- function(spec, axis) {
  if (!axis) return(0L)
  if (is.null(spec$gradient.order)) 1L else spec$gradient.order[axis]
}

# d^r exp(-(x-t)^2/(2h^2))/dt^r divided by that Gaussian weight.
.np_bias_gaussian_derivatives <- function(u, h, order) {
  out <- matrix(1,length(u),order+1L)
  if (order>0L) out[,2L] <- u/h
  if (order>1L) for (r in 2:order)
    out[,r+1L] <- (u*out[,r]-((r-1)/h)*out[,r-1L])/h
  out
}

.np_bias_weighted_derivative <- function(Q, weight, value, ratios,
                                         penalty=NULL, role="higher-order pilot", rank.tolerance=1e-12) {
  P <- ncol(Q); n <- nrow(Q); order <- nrow(value)-1L
  if (!is.matrix(value) || ncol(value)!=P || order<1L ||
      !is.matrix(ratios) || !identical(dim(ratios),c(n,order+1L)) ||
      length(weight)!=n || any(!is.finite(weight)) || any(weight<0))
    stop("invalid higher-order weighted-pilot geometry",call.=FALSE)
  if (is.null(penalty)) penalty <- numeric(P)
  if (length(penalty)!=P || any(!is.finite(penalty)) || any(penalty<0))
    stop("invalid fixed auxiliary penalty",call.=FALSE)
  weighted <- Q*sqrt(weight)
  augmented <- if (any(penalty>0)) rbind(weighted,diag(sqrt(penalty),P)) else weighted
  norm <- sqrt(colSums(augmented^2))
  if (any(!is.finite(norm)) || any(norm<=0)) stop(role," has an unsupported column",call.=FALSE)
  fit <- qr(sweep(augmented,2L,norm,"/"),tol=rank.tolerance,LAPACK=FALSE)
  if (fit$rank!=P) stop(role," is numerically unidentified",call.=FALSE)
  R <- qr.R(fit); pivot <- fit$pivot
  solve.gram <- function(rhs) {
    z <- backsolve(R,forwardsolve(t(R),rhs[pivot,,drop=FALSE]))
    out <- z; out[pivot,] <- z
    out
  }
  Z <- sweep(Q,2L,norm,"/"); V <- sweep(value,2L,norm,"/")
  A <- G <- vector("list",order+1L)
  A[[1L]] <- solve.gram(t(Z)*rep(weight,each=P))
  for (r in seq_len(order)) {
    w <- weight*ratios[,r+1L]
    G[[r+1L]] <- crossprod(Z,Z*w)
    rhs <- t(Z)*rep(w,each=P)
    for (j in seq_len(r)) rhs <- rhs-choose(r,j)*(G[[j+1L]]%*%A[[r-j+1L]])
    A[[r+1L]] <- solve.gram(rhs)
  }
  out <- matrix(0,1L,n)
  for (j in 0:order) out <- out+choose(order,j)*(V[j+1L,,drop=FALSE]%*%A[[order-j+1L]])
  if (any(!is.finite(out))) stop(role," produced nonfinite derivative rows",call.=FALSE)
  drop(out)
}
