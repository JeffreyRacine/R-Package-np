# Qualified local fixed-bandwidth calibration. Other inference owners retain
# their existing variance and main-row contracts.
.np_bias_local_supported <- function(state, law) {
  if (!is.character(law) || length(law)!=1L ||
      !law %in% c("gaussian","mammen") || state$d!=1L ||
      length(state$factor.cols) || isTRUE(state$pooled) ||
      isTRUE(state$categorical.only) || !is.null(state$adaptation) ||
      !is.null(state$mixed) ||
      !.np_plot_regression_capture_supported(state$bws) ||
      (!is.null(state$gradient.order) && any(state$gradient.order!=1L)))
    return(FALSE)
  reg <- npValidatedConditionalRegSpec(state$bws,
    where="local bias calibration",ncon.field="ncon")
  identical(reg[["basis.engine"]],"glp") &&
    (all(reg[["degree.engine"]]<=1L) || isTRUE(reg[["bernstein.basis.engine"]]))
}

# Design-only neighborhood geometry, storing O(n) radii and normalizers.
.np_bias_local_neighborhood <- function(x, q) {
  x <- as.matrix(x); n <- nrow(x); d <- ncol(x)
  scales <- apply(x,2L,sd)
  if (n<2L || d<1L || length(q)!=n || any(!is.finite(q)) || any(q<=0) ||
      any(!is.finite(x)) || any(!is.finite(scales)) || any(scales<=0))
    stop("invalid local calibration neighborhood geometry",call.=FALSE)
  z <- sweep(sweep(x,2L,colMeans(x),"-"),2L,scales,"/")
  k <- as.integer(ceiling(n^(4/(4+d))))
  distance <- function(ii) {
    if (!length(ii) || length(ii)>16L || anyNA(ii) ||
        any(ii<1L | ii>n | ii!=floor(ii)))
      stop("invalid local calibration neighborhood tile",call.=FALSE)
    D <- matrix(0,length(ii),n)
    for (j in seq_len(d)) D <- D+outer(z[ii,j],z[,j],"-")^2
    D
  }
  radius2 <- denominator <- count <- numeric(n)
  for (first in seq.int(1L,n,by=16L)) {
    ii <- seq.int(first,min(n,first+15L)); D <- distance(ii)
    radius2[ii] <- apply(D,1L,function(v)sort(v,partial=k)[k])
    selected <- D<=radius2[ii]+64*.Machine$double.eps*pmax(1,radius2[ii])
    count[ii] <- rowSums(selected); denominator[ii] <- drop(selected%*%q)
  }
  if (any(count<k) || any(!is.finite(denominator)) || any(denominator<=0))
    stop("invalid local calibration neighborhood mass",call.=FALSE)
  raw <- function(ii)
    (distance(ii)<=radius2[ii]+64*.Machine$double.eps*pmax(1,radius2[ii]))/denominator[ii]
  list(raw=raw,k=k,radius2=radius2,count=count,scales=scales)
}

.np_bias_local_prepare <- function(state) {
  q <- state$residual_row_variance; raw <- state$residual_raw
  n <- length(q)
  smoother <- list(x=state$xcon,groups=state$groups,h=state$variance_h)
  S <- .np_bias_student_variance_owner(smoother)
  neighborhood <- .np_bias_local_neighborhood(state$xcon,q)
  denominator <- variance <- generation <- numeric(n)
  for (first in seq.int(1L,n,by=16L)) {
    ii <- seq.int(first,min(n,first+15L)); H <- S(ii)
    denominator[ii] <- drop(H%*%q)
    variance[ii] <- drop(H%*%raw^2)/denominator[ii]
    generation[ii] <- drop(neighborhood$raw(ii)%*%raw^2)
    .np_plot_progress_notify(done=max(ii),total=n,unit="training rows")
  }
  if (any(!is.finite(c(denominator,variance,generation))) ||
      any(denominator<=0) || any(variance<=0) || any(generation<=0))
    stop("nonpositive local bias calibration variance",call.=FALSE)
  state$variance <- variance
  state$local.inference <- list(version="captured-pooled-positive-v1",
    denominator=denominator,generation=generation)
  state$settings$local.inference <- list(version="captured-pooled-positive-v1",
    main="complete original fitted capture",observed="S(raw residual^2) / S(row variance)",
    generation="positive normalized-residual neighborhood variance",
    k=neighborhood$k,count=neighborhood$count,radius2=neighborhood$radius2,
    coordinate.scale=neighborhood$scales,
    qualification="univariate continuous fixed Gaussian2, inactive adaptation, mean/first derivative, Gaussian/Mammen")
  state
}
