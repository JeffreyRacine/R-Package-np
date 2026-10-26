# Private fixed-density variance pilot. All moments use the same h/lambda.
# The estimating center and existing empirical-covariance adapter are untouched.
.np_udens_b2_geosums <- function(r, n = Inf) {
  if (is.finite(n) && n <= 64) {
    j <- seq_len(n) - 1
    return(c(sum(r^j), sum(j*r^j), sum(j^2*r^j)))
  }
  if (r == 1) return(c(n, n*(n-1)/2, n*(n-1)*(2*n-1)/6))
  a <- c(1/(1-r), r/(1-r)^2, r*(1+r)/(1-r)^3)
  if (!is.finite(n)) return(a)
  a - r^n*c(a[1], n*a[1]+a[2], n^2*a[1]+2*n*a[2]+a[3])
}

.np_udens_b2_lattice_moment <- function(donor, target, lambda, corrected) {
  # The pilot S is normalized on Z; target contains one or two Q arguments.
  c <- (1-lambda)/(1+lambda)
  A <- if (corrected) 2-c-2*lambda^2/(1+lambda)^2 else 1
  slope <- if (corrected) c else 0
  if (lambda == 0) return(as.double(all(target == donor)))
  knots <- sort(unique(c(donor, target)))
  lo <- c(-Inf, knots+1); hi <- c(knots, Inf)
  total <- 0
  for (j in seq_along(lo)) {
    if (lo[j] > hi[j]) next
    exponent <- function(z) abs(z-donor)+sum(abs(z-target))
    if (!is.finite(lo[j])) {z <- hi[j]; direction <- -1; count <- Inf} else
    if (!is.finite(hi[j])) {z <- lo[j]; direction <- 1; count <- Inf} else {
      count <- hi[j]-lo[j]+1
      direction <- if (exponent(lo[j]) <= exponent(hi[j])) 1 else -1
      z <- if (direction == 1) lo[j] else hi[j]
    }
    e <- exponent(z)
    rate <- if (count == 1) 0 else exponent(z+direction)-e
    if (rate < 0) stop("invalid lattice moment partition", call. = FALSE)
    u <- A-slope*abs(z-target)
    v <- -slope*(abs(z+direction-target)-abs(z-target))
    coefficient <- if (length(target)==1L) c(u,v,0) else
      c(u[1]*u[2],u[1]*v[2]+u[2]*v[1],v[1]*v[2])
    total <- total + c^(length(target)+1L)*lambda^e*
      sum(coefficient*.np_udens_b2_geosums(lambda^rate,count))
  }
  total
}

.np_udens_b2_category_moment <- function(category, target, corrected = TRUE) {
  donors <- unique(category$training)
  values <- vapply(donors, function(donor) {
    if (category$ordered)
      return(.np_udens_b2_lattice_moment(donor,target,category$lambda,corrected))
    k <- category$k; lambda <- category$lambda
    labels <- unique(c(donor,target)); extra <- k-length(labels)
    if (extra < 0) stop("categorical moment support mismatch",call.=FALSE)
    if (k == 1L) return(1)
    pilot <- ifelse(labels==donor,1-lambda,lambda/(k-1))
    eff <- if (corrected) k*lambda^2/(k-1) else lambda
    product <- rep(1,length(labels))
    for (t in target) product <- product*ifelse(labels==t,1-eff,eff/(k-1))
    sum(pilot*product)+extra*lambda/(k-1)*(eff/(k-1))^length(target)
  },numeric(1L))
  values[match(category$training,donors)]
}

.np_udens_b2_pilot_moment <- function(state, left, right = NULL, corrected = TRUE) {
  numeric <- vapply(state$x,is.numeric,logical(1L));q <- sum(numeric)
  x <- state$x[numeric];h <- state$h[numeric];t <- as.double(left[numeric])
  weight <- rep(1,state$n);et <- es <- cross <- numeric(state$n)
  if (is.null(right)) {
    radius <- numeric(state$n)
    for (j in seq_len(q)) {
      u <- (t[j]-x[[j]])/h[j]
      weight <- weight*dnorm(u/sqrt(2))/(sqrt(2)*h[j]);radius <- radius+u^2
    }
    value <- weight*if(corrected) (1+q/4-radius/8) else 1
  } else {
    s <- as.double(right[numeric])
    for (j in seq_len(q)) {
      a <- (t[j]-x[[j]])/h[j];b <- (s[j]-x[[j]])/h[j]
      weight <- weight*exp(-((a-b)^2+a^2+b^2)/6)/(2*pi*h[j]^2*sqrt(3))
      mt <- (2*a-b)/3;ms <- (2*b-a)/3
      et <- et+mt^2+1/3;es <- es+ms^2+1/3;cross <- cross+4*mt*ms/3+2/9
    }
    value <- weight*if(corrected) ((1+q/2-et/2)*(1+q/2-es/2)+cross/4) else 1
  }
  if (state$r+state$s > 0) for (z in state$category) {
    target <- as.double(toMatrix(left[z$index]))
    if (!is.null(right)) target <- c(target,as.double(toMatrix(right[z$index])))
    value <- value*.np_udens_b2_category_moment(z,target,corrected)
  }
  if (any(!is.finite(value))) stop("nonfinite smoothed density moment",call.=FALSE)
  value
}

.np_udens_b2_smoothed_covariance <- function(state, grid, tile = 25L, corrected = TRUE) {
  .np_udens_bc_grid(state,grid)
  m <- nrow(grid);.np_udens_bc_budget(state,m,1L,tile,TRUE)
  mu <- center <- original <- numeric(m);C <- matrix(0,m,m)
  for (block in .np_bias_target_blocks(m,1L,anchors=tile)) {
    ii <- block$grid;z <- .np_udens_bc_contribution(state,grid[ii,,drop=FALSE],corrected)
    center[ii] <- z$center;original[ii] <- z$original.center
  }
  for (i in seq_len(m)) {
    mu[i] <- mean(.np_udens_b2_pilot_moment(state,grid[i,,drop=FALSE],corrected=corrected))
    for (j in seq_len(i)) C[i,j] <- C[j,i] <-
      mean(.np_udens_b2_pilot_moment(state,grid[i,,drop=FALSE],grid[j,,drop=FALSE],corrected))
  }
  C <- (C-tcrossprod(mu))/state$n
  if (any(!is.finite(C)) || any(diag(C)<=0))
    stop("smoothed density covariance requires positive finite scales",call.=FALSE)
  list(center=center,original.center=original,se=sqrt(diag(C)),covariance=C,
       pilot.mean=mu,variance.measure="same-bandwidth normalized positive KDE")
}

# Reuse the existing bounded Gaussian engine by replacing only its private
# moment dependency in local closures. No namespace binding is changed.
.np_udens_b2_gaussian <- function(state, grid, alpha = c(.05,.1), B = 9999L,
                                  tile = 25L, corrected = TRUE) {
  owner <- new.env(parent = environment(.np_udens_bc_gaussian))
  owner$.np_udens_bc_covariance <- .np_udens_b2_smoothed_covariance
  engine <- .np_udens_bc_gaussian; environment(engine) <- owner
  owner$.np_udens_bc_gaussian <- engine
  entry <- .np_udens_bc_inference; environment(entry) <- owner
  out <- entry(state,grid,law="gaussian",alpha=alpha,B=B,tile=tile,corrected=corrected)
  out$settings$construction <- "density-bc-v2-smoothed-covariance"
  out$settings$covariance <- "same-bandwidth positive-KDE measure, n-squared convention"
  out
}
