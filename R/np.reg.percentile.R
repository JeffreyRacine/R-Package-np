# Experimental E1-R mean / E6-R derivative calibrated percentiles.
# Original fitted geometry is immutable. Empirical standardized innovations
# are resampled at fixed design; no kernel jitter is added to residuals.
.np_reg_percentile_with_seed <- function(seed, code) {
 had <- exists('.Random.seed',.GlobalEnv,inherits=FALSE)
 old <- if(had)get('.Random.seed',.GlobalEnv) else NULL
 kind <- RNGkind()
 on.exit({do.call(RNGkind,as.list(kind));if(had)assign('.Random.seed',old,.GlobalEnv) else if(exists('.Random.seed',.GlobalEnv,inherits=FALSE))rm('.Random.seed',envir=.GlobalEnv)})
 set.seed(seed);force(code)
}
# Capture opaque runtime state without selecting a different generator.
.np_reg_percentile_rng_mode <- function() {
  if (!identical(RNGkind()[1L], 'user-supplied')) return('state')
  .npRmpi_bootstrap_fanout_enabled(what = 'regression percentile calibration')
  if (isTRUE(.npRmpi_has_active_slave_pool()) ||
      !isTRUE(.npRmpi_master_only_mode()))
    stop('regression smoothed-percentile cannot reproduce a user-supplied uniform RNG on workers; use a built-in RNG or an explicitly selected supported local execution route', call. = FALSE)
  'legacy-local-custom'
}
.np_reg_percentile_rng_state <- function(seed) {
  had <- exists('.Random.seed', .GlobalEnv, inherits = FALSE)
  old <- if (had) get('.Random.seed', .GlobalEnv, inherits = FALSE) else NULL
  on.exit(if (had) assign('.Random.seed', old, .GlobalEnv) else
    if (exists('.Random.seed', .GlobalEnv, inherits = FALSE))
      rm('.Random.seed', envir = .GlobalEnv))
  set.seed(seed)
  get('.Random.seed', .GlobalEnv, inherits = FALSE)
}
.np_reg_percentile_rng_validate <- function(state, contract) {
  invalid <- function() stop('invalid regression percentile task RNG state', call. = FALSE)
  kind.length <- if (getRversion() < '3.6.0') 2L else 3L
  if (!is.integer(state) || length(state) < 2L || is.na(state[1L]) ||
      !is.list(contract) ||
      !identical(contract$runtime, paste(R.version$major, R.version$minor, sep = '.')) ||
      !identical(contract$header, state[1L]) ||
      !identical(contract$length, length(state)) ||
      !is.character(contract$kind) || length(contract$kind) != kind.length ||
      anyNA(contract$kind)) invalid()
  # R_ext/Random.h kind codes; state tail words remain opaque integer bits.
  uniform <- match(contract$kind[1L], c('Wichmann-Hill', 'Marsaglia-Multicarry',
    'Super-Duper', 'Mersenne-Twister', 'Knuth-TAOCP', 'user-supplied',
    'Knuth-TAOCP-2002', "L'Ecuyer-CMRG")) - 1L
  normal <- match(contract$kind[2L], c('Buggy Kinderman-Ramage', 'Ahrens-Dieter',
    'Box-Muller', 'user-supplied', 'Inversion', 'Kinderman-Ramage')) - 1L
  sample <- if (kind.length == 2L) 0L else
    match(contract$kind[3L], c('Rounding', 'Rejection')) - 1L
  if (anyNA(c(uniform, normal, sample)) || uniform == 5L ||
      !identical(state[1L], uniform + 100L * normal + 10000L * sample) ||
      length(state) != c(4L, 3L, 3L, 626L, 102L, 0L, 102L, 7L)[uniform + 1L]) invalid()
  invisible(TRUE)
}
.np_reg_percentile_with_state <- function(state, contract, code) {
  .np_reg_percentile_rng_validate(state, contract)
  had <- exists('.Random.seed', .GlobalEnv, inherits = FALSE)
  old <- if (had) get('.Random.seed', .GlobalEnv, inherits = FALSE) else NULL
  kind <- RNGkind(); complete <- FALSE
  on.exit({
    warnings <- list()
    # Restore before notifying caller handlers, which may themselves throw.
    cleanup <- tryCatch(withCallingHandlers(do.call(RNGkind, as.list(kind)),
      warning = function(w) {
        warnings[[length(warnings) + 1L]] <<- w
        invokeRestart('muffleWarning')
      }), error = identity, interrupt = identity)
    if (had) assign('.Random.seed', old, .GlobalEnv) else
      if (exists('.Random.seed', .GlobalEnv, inherits = FALSE))
        rm('.Random.seed', envir = .GlobalEnv)
    # Never replace a condition already unwinding through the draw body.
    if (complete) {
      if (inherits(cleanup, 'condition')) stop(cleanup)
      for (w in warnings) warning(w)
    }
  })
  assign('.Random.seed', state, .GlobalEnv)
  value <- force(code)
  complete <- TRUE
  value
}
.np_reg_percentile_rank_bands <- function(curves,map,alpha=.05) {
 B<-ncol(curves); families<-unique(map$target)
 bands<-setNames(lapply(1:3,function(i)array(0,c(nrow(map),2L,length(alpha)))),c('pointwise','simultaneous','bonferroni'))
 details<-list()
 for(g in families) {
  ii<-which(map$target==g); mm<-map[ii,,drop=FALSE]
  z<-.np_band_rank_prepare(B,mm,alpha)
  z<-.np_band_rank_update(z,curves[ii,,drop=FALSE],mm,seq_len(B))
  z<-.np_band_rank_cutoff(z)
  z<-.np_band_rank_update(z,curves[ii,,drop=FALSE],mm,seq_len(B),replay=TRUE)
  z<-.np_band_rank_finish(z)
  z$orientation<-'direct corrected response curves, no inversion or recentering'
  for(b in names(bands))bands[[b]][ii,,]<-z$limits[[b]]
  details[[g]]<-z[c('kP','kS','kB','content','ties','orientation')]
 }
 # Resolution helper orders families lexicographically.
 ks<-do.call(rbind,lapply(sort(families),function(g)details[[g]]$kS))
 resolution<-.np_band_resolution(B,map$target,alpha,rep(TRUE,nrow(map)),names(bands),signed=TRUE,simultaneous.index=ks)
 resolution$calibration<-'smoothed-percentile'
 .np_band_resolution_warn(resolution)
 list(bands=bands,ranks=details,resolution=resolution)
}
.np_reg_percentile_shape <- function(p) {
 v<-p$residual/(sqrt(p$row.variance)*p$sigma)
 v<-v-mean(v); den<-sqrt(mean(v^2))
 if(any(!is.finite(v))||!is.finite(den)||den<=0)stop('undefined empirical residual normalization')
 u<-v/den
 list(u=u,third=mean(u^3),fourth=mean(u^4),mean=mean(u),variance=mean(u^2))
}
.np_reg_percentile_refresh_geometry <- function(p) {
 sp<-.np_bias_student_prepare(p$state)
 H<-.np_bias_student_pilot_owner(sp$pilot)
 local<-!is.null(p$state$local.inference)
 if(local) {
  ng<-.np_bias_local_neighborhood(p$state$xcon,p$row.variance)
  variance<-ng$raw
 } else variance<-.np_bias_student_variance_owner(sp$smoother)
 list(H=H,variance=variance,local=local,q=p$row.variance)
}
.np_reg_percentile_rank_info <- function(x){
 stopifnot(is.matrix(x),all(is.finite(x)),ncol(x)>=1,nrow(x)>=1)
 B<-ncol(x);sorted<-matrix(0,nrow(x),B);K<-rep.int(B,B)
 for(i in seq_len(nrow(x))){
  o<-order(x[i,],method='radix');v<-x[i,o];sorted[i,]<-v
  runs<-rle(v)$lengths;ends<-cumsum(runs);starts<-ends-runs+1L
  rk<-integer(B);rk[o]<-rep(pmin(ends,B+1L-starts),runs);K<-pmin(K,rk)
 }
 list(sorted=sorted,K=K,ordered=sort(K,method='radix'),B=B,m=nrow(x))
}
.np_reg_percentile_companion_k <- function(info,j){
 B<-info$B;stopifnot(length(j)==1,j>=1,j<=B,j==as.integer(j))
 min(1+floor(B*.05/2),1+floor((B-j)/2),info$ordered[B-j+1L])
}
.np_reg_percentile_least_content <- function(info,truth){
 stopifnot(length(truth)==info$m,all(is.finite(truth)))
 d<-min(vapply(seq_len(info$m),function(i)min(sum(info$sorted[i,]<=truth[i]),sum(info$sorted[i,]>=truth[i])),0))
 B<-info$B;j0<-ceiling(.501*B);kp<-1+floor(B*.05/2)
 if(d>=kp)return(j0)
 if(d==0)return(Inf)
 j<-max(j0,min(B-2*d+1,B-sum(info$K<=d)+1))
 if(j>B)return(Inf)
 stopifnot(.np_reg_percentile_companion_k(info,j)<=d)
 j
}
.np_reg_percentile_companion_bands <- function(info,j,nominal){
 m<-info$m;B<-info$B;k<-if(is.finite(j)).np_reg_percentile_companion_k(info,j) else 1L
 S<-array(c(info$sorted[,k],info$sorted[,B+1-k]),c(m,2L,1L))
 out<-nominal
 out$simultaneous<-S
 out$simultaneous[,1,1]<-pmin(nominal$pointwise[,1,1],S[,1,1])
 out$simultaneous[,2,1]<-pmax(nominal$pointwise[,2,1],S[,2,1])
 out$bonferroni[,1,1]<-pmin(nominal$bonferroni[,1,1],out$simultaneous[,1,1])
 out$bonferroni[,2,1]<-pmax(nominal$bonferroni[,2,1],out$simultaneous[,2,1])
 if(m==1L)out$pointwise<-out$bonferroni<-out$simultaneous
 out
}
.np_reg_percentile_least_content_fast<-function(x,truth){
 stopifnot(is.matrix(x),nrow(x)>=1L,ncol(x)>=1L,all(is.finite(x)),length(truth)==nrow(x),all(is.finite(truth)))
 B<-ncol(x);d<-min(vapply(seq_len(nrow(x)),function(i)min(sum(x[i,]<=truth[i]),sum(x[i,]>=truth[i])),0))
 j0<-ceiling(.501*B);kp<-1+floor(B*.05/2)
 if(d>=kp)return(as.double(j0))
 if(d==0)return(Inf)
 outside<-rep(FALSE,B)
 for(i in seq_len(nrow(x))){
  v<-x[i,];vsel<-sort(v,partial=c(d+1L,B-d))
  outside<-outside | v<vsel[d+1L] | v>vsel[B-d]
 }
 as.double(max(j0,min(B-2*d+1,B-sum(outside)+1)))
}
.np_reg_percentile_e6_scores <- function(x, truth) {
 stopifnot(is.matrix(x),nrow(x)>0L,ncol(x)>0L,all(is.finite(x)),
           length(truth)==nrow(x),all(is.finite(truth)))
 B<-ncol(x);k0<-ceiling(B/2)
 depths<-vapply(seq_len(nrow(x)),function(i)min(sum(x[i,]<=truth[i]),sum(x[i,]>=truth[i])),0)
 marginal<-ifelse(depths>=k0,0,ifelse(depths==0,Inf,pmax(1,B-2*depths+1)))
 d<-min(depths)
 family<-if(d>=k0)0 else if(d==0)Inf else {
  outside<-rep(FALSE,B)
  for(i in seq_len(nrow(x))){v<-x[i,];s<-sort(v,partial=c(d+1L,B-d));outside<-outside|v<s[d+1L]|v>s[B-d]}
  max(1,min(B-2*d+1,B-sum(outside)+1))
 }
 list(marginal=as.double(marginal),family=as.double(family))
}
.np_reg_percentile_e6_index <- function(B,j) {
 stopifnot(all(!is.na(j)),all(j>=0),all(!is.finite(j)|j<=B),all(!is.finite(j)|j==floor(j)))
 ifelse(!is.finite(j),1,ifelse(j==0,ceiling(B/2),1+floor((B-j)/2)))
}
.np_reg_percentile_e6_family_index <- function(info,j) {
 if(!is.finite(j))return(1)
 if(j==0)return(ceiling(info$B/2))
 min(.np_reg_percentile_e6_index(info$B,j),info$ordered[info$B-j+1L])
}
.np_reg_percentile_e6_materialize <- function(info,marginal,family,support,nominal_bonf) {
 B<-info$B;m<-info$m;stopifnot(length(marginal)==m,all(dim(support)==c(m,2L)))
 kp<-.np_reg_percentile_e6_index(B,marginal);ks<-.np_reg_percentile_e6_family_index(info,family)
 P<-cbind(info$sorted[cbind(seq_len(m),kp)],info$sorted[cbind(seq_len(m),B+1-kp)])
 S<-cbind(info$sorted[,ks],info$sorted[,B+1-ks])
 boundary_p<-!is.finite(marginal)|marginal==B;boundary_s<-!is.finite(family)||family==B
 P[boundary_p,1]<-pmin(P[boundary_p,1],support[boundary_p,1]);P[boundary_p,2]<-pmax(P[boundary_p,2],support[boundary_p,2])
 if(boundary_s){S[,1]<-pmin(S[,1],support[,1]);S[,2]<-pmax(S[,2],support[,2])}
 S[,1]<-pmin(S[,1],P[,1]);S[,2]<-pmax(S[,2],P[,2])
 C<-cbind(pmin(nominal_bonf[,1],S[,1]),pmax(nominal_bonf[,2],S[,2]))
 if(m==1L)P<-C<-S
 list(bands=list(pointwise=array(P,c(m,2L,1L)),simultaneous=array(S,c(m,2L,1L)),
                 bonferroni=array(C,c(m,2L,1L))),kp=kp,ks=ks,
      support.used=any(boundary_p)||boundary_s,
      attainable=all(is.finite(marginal))&&is.finite(family))
}
.np_reg_percentile_outer_response_bank<-function(p,g,outer_data,law){
 stopifnot(length(outer_data)==39L)
 ys<-lapply(outer_data,function(seed).np_reg_percentile_with_seed(seed,p$mu+p$sigma*drop(.np_reg_percentile_draw_innovations(p,1L))))
 n<-length(p$y);raw<-v<-matrix(0,n,39L)
 for(first in seq.int(1L,n,16L)){
  ii<-seq.int(first,min(n,first+15L));H<-g$H(ii)
  for(o in 1:39)raw[ii,o]<-ys[[o]][ii]-drop(H%*%ys[[o]])
 }
 squares<-raw
 for(o in 1:39)squares[,o]<-if(g$local)raw[,o]^2 else raw[,o]^2/g$q
 for(first in seq.int(1L,n,16L)){
  ii<-seq.int(first,min(n,first+15L));V<-g$variance(ii)
  for(o in 1:39)v[ii,o]<-drop(V%*%squares[,o])
 }
 for(o in 1:39)if(any(!is.finite(v[,o]))||any(v[,o]<=0))stop('undefined outer generation variance')
 list(y=ys,raw=raw,v=v)
}
.np_reg_percentile_bank_refreshed<-function(p,y,bank,o){
 q<-p;q$y<-y;q$mu<-y-bank$raw[,o];q$sigma<-sqrt(bank$v[,o]);q$residual<-bank$raw[,o]
 q$center<-drop(p$A%*%y);q
}

.np_reg_percentile_draw_innovations <- function(p, size) {
  u <- .np_reg_percentile_shape(p)$u; n <- length(u)
  matrix(u[floor(n * runif(n * size)) + 1L], n, size)
}
.np_reg_percentile_curve_draws <- function(p, B, seed, rng.mode = 'seed',
                                             rng.state = NULL, rng.contract = NULL) {
  draw <- function() {
  curves <- matrix(0, nrow(p$A), B)
  for (first in seq.int(1L, B, 256L)) {
    jj <- seq.int(first, min(B, first + 255L))
    ep <- p$sigma * .np_reg_percentile_draw_innovations(p, length(jj))
    curves[, jj] <- p$A %*% sweep(ep, 1L, p$mu, '+')
  }
  list(curves = curves, end = .Random.seed)
  }
  if (identical(rng.mode, 'state'))
    return(.np_reg_percentile_with_state(rng.state, rng.contract, draw()))
  if (!rng.mode %in% c('seed', 'legacy-local-custom'))
    stop('invalid regression percentile task RNG mode', call. = FALSE)
  .np_reg_percentile_with_seed(seed, draw())
}
.np_reg_percentile_curve_support <- function(p) {
  u <- range(.np_reg_percentile_shape(p)$u); a <- sweep(p$A, 2L, p$sigma, '*')
  mu <- drop(p$A %*% p$mu)
  L <- mu + rowSums(pmin(a * u[1L], a * u[2L]))
  U <- mu + rowSums(pmax(a * u[1L], a * u[2L]))
  if (any(!is.finite(c(L, U))) || any(L > U))
    stop('undefined regression percentile support', call. = FALSE)
  cbind(L, U)
}
.np_reg_percentile_outer_score <- function(task, p, B) {
  p$mu <- task$mu; p$sigma <- task$sigma; p$residual <- task$raw
  if (!identical(task$rng.mode, 'state') &&
      !identical(task$rng.mode, 'legacy-local-custom'))
    stop('missing regression percentile task RNG mode', call. = FALSE)
  curves <- .np_reg_percentile_curve_draws(p, B, task$seed,
    task$rng.mode, task$rng.state, task$rng.contract)$curves
  support <- .np_reg_percentile_curve_support(p)
  scores <- empirical <- numeric()
  for (target in unique(p$map$target)) {
    ii <- which(p$map$target == target); truth <- p$pilot.eval[ii]
    X <- curves[ii, , drop = FALSE]
    contained <- support[ii, 1L] <= truth & truth <= support[ii, 2L]
    if (target == 'mean') {
      value <- .np_reg_percentile_least_content_fast(X, truth); complete <- all(contained)
    } else {
      z <- .np_reg_percentile_e6_scores(X, truth); value <- c(z$marginal, z$family)
      complete <- c(contained, all(contained))
    }
    empirical <- c(empirical, value)
    value[!is.finite(value) & complete] <- B
    scores <- c(scores, value)
  }
  # Finite transport sentinel only: decode before any statistical operation.
  z <- c(scores, empirical); z[is.infinite(z)] <- B + 1
  matrix(z, 1L)
}
# Coordinator prepares geometry once; workers own disjoint complete outer curves.
.np_reg_percentile_outer_dispatch <- function(tasks, p, B) {
  .npRmpi_bootstrap_fanout_enabled(n=ncol(p$A), B=39L, chunk.size=1L,
    what="regression percentile calibration")
  .np_reg_percentile_rng_mode()
  if (isTRUE(.npRmpi_has_active_slave_pool()) &&
      !all(vapply(tasks, function(task) identical(task$rng.mode, 'state'), logical(1L))))
    stop('distributed regression percentile requires explicit task RNG states', call. = FALSE)
  for (task in tasks) if (identical(task$rng.mode, 'state'))
    .np_reg_percentile_rng_validate(task$rng.state, task$rng.contract)
  bindings <- mget(c(".np_reg_percentile_with_state", ".np_reg_percentile_rng_validate",
    ".np_reg_percentile_outer_score", ".np_reg_percentile_curve_draws",
    ".np_reg_percentile_curve_support", ".np_reg_percentile_with_seed",
    ".np_reg_percentile_draw_innovations", ".np_reg_percentile_shape",
    ".np_reg_percentile_least_content_fast", ".np_reg_percentile_e6_scores"),
    envir=environment(.np_reg_percentile_outer_dispatch), inherits=TRUE)
  size <- 1L + sum(vapply(setdiff(unique(p$map$target),"mean"),
    function(target) sum(p$map$target==target)+1L, integer(1L)))
  worker <- function(task, percentile.state, B)
    .np_reg_percentile_outer_score(task, percentile.state, B)
  .npRmpi_bootstrap_run_fanout(tasks, worker,
    ncol.out=2L*size, what="regression percentile calibration",
    progress.label="Calibrating regression percentile outer samples",
    profile.where="mpi.applyLB:regression-percentile", required.bindings=bindings,
    percentile.state=p, B=B, progress.unit="samples")
}
.np_reg_percentile_forecast <- function(n, m, B) {
  if (any(!is.finite(c(n,m,B))) || any(c(n,m,B)<1) ||
      any(c(n,m,B)!=floor(c(n,m,B))))
    stop("invalid regression percentile forecast dimensions", call.=FALSE)
  n <- as.double(n); m <- as.double(m); B <- as.double(B)
  bytes <- 8 * (10 * m * B + 10 * n * 39 + 6 * m * n + 20 * n * 256)
  if (!is.finite(bytes) || bytes > 256 * 1024^2)
    stop("regression percentile memory forecast exceeds 256 MiB; reduce the requested grid or B", call.=FALSE)
  bytes
}

.np_reg_percentile_infer <- function(p, B, seed, outer_data, outer_inner) {
  B <- .np_plot_gaussian_integer(B, 'regression percentile draws', 9999L)
  stopifnot(length(outer_data) == 39L, length(outer_inner) == 39L)
  n <- length(p$y); m <- nrow(p$A)
  # Original and current banks, sorting temporaries, bounded geometry and tasks.
  forecast <- .np_reg_percentile_forecast(n, m, B)
  rng.mode <- .np_reg_percentile_rng_mode()
  g <- .np_reg_percentile_refresh_geometry(p)
  bank <- .np_reg_percentile_outer_response_bank(p, g, outer_data, 'residual')
  original <- .np_reg_percentile_curve_draws(p, B, seed)$curves
  nominal <- suppressWarnings(.np_reg_percentile_rank_bands(original, p$map))
  tasks <- lapply(seq_len(39L), function(o) list(start = o, bsz = 1L,
    seed = outer_inner[o], mu = bank$y[[o]] - bank$raw[, o],
    sigma = sqrt(bank$v[, o]), raw = bank$raw[, o]))
  for (o in seq_along(tasks)) {
    tasks[[o]]$rng.mode <- rng.mode
    if (identical(rng.mode, 'state')) {
      state <- .np_reg_percentile_rng_state(outer_inner[o])
      tasks[[o]]$rng.state <- state
      tasks[[o]]$rng.contract <- list(
        runtime = paste(R.version$major, R.version$minor, sep = '.'),
        header = state[1L], length = length(state), kind = RNGkind())
    }
  }
  payload <- p[c('A', 'mu', 'sigma', 'residual', 'row.variance', 'pilot.eval', 'map')]
  rows <- .np_plot_progress_outer_run(
    .np_reg_percentile_outer_dispatch(tasks, payload, B))
  rows[rows == B + 1] <- Inf
  size <- ncol(rows) / 2L
  scores <- rows[, seq_len(size), drop = FALSE]
  empirical <- rows[, size + seq_len(size), drop = FALSE]
  support <- .np_reg_percentile_curve_support(p); z <- nominal; z$calibration <- list(); at <- 0L
  for (target in unique(p$map$target)) {
    ii <- which(p$map$target == target); m <- length(ii)
    span <- if (target == 'mean') 1L else m + 1L
    jj <- at + seq_len(span); at <- at + span
    ss <- scores[, jj, drop = FALSE]; ee <- empirical[, jj, drop = FALSE]
    sel <- apply(ss, 2L, function(v) sort(v)[38L])
    info <- .np_reg_percentile_rank_info(original[ii, , drop = FALSE])
    if (target == 'mean') {
      s <- sel[1L]; small <- lapply(nominal$bands, function(x) x[ii, , , drop = FALSE])
      bands <- .np_reg_percentile_companion_bands(info, s, small); boundary <- !is.finite(s) || s == B
      if (boundary) {
        L <- pmin(support[ii, 1L], info$sorted[, 1L], bands$simultaneous[, 1L, 1L])
        U <- pmax(support[ii, 2L], info$sorted[, B], bands$simultaneous[, 2L, 1L])
        bands$simultaneous[, 1L, 1L] <- L; bands$simultaneous[, 2L, 1L] <- U
        bands$bonferroni[, 1L, 1L] <- pmin(bands$bonferroni[, 1L, 1L], L)
        bands$bonferroni[, 2L, 1L] <- pmax(bands$bonferroni[, 2L, 1L], U)
        if (m == 1L) bands$pointwise <- bands$bonferroni <- bands$simultaneous
      }
      z$calibration[[target]] <- list(scores = ss[, 1L], empirical.scores = ee[, 1L],
        selected = s, q = s / B, rank = 38L, J = 39L, attainable = is.finite(s),
        support.used = boundary, outer.completions = sum(!is.finite(ee) & is.finite(ss)),
        floor = sum(ss == ceiling(.501 * B)), support = support[ii, , drop = FALSE], law = 'residual')
      ks <- if (is.finite(s)) .np_reg_percentile_companion_k(info, s) else 1
      kp <- nominal$ranks[[target]]$kP
      kb <- min(nominal$ranks[[target]]$kB, ks)
      if (m == 1L) kp <- kb <- ks
      z$ranks[[target]][c('kP', 'kS', 'kB', 'content')] <-
        list(kp, ks, kb, sum(info$K >= ks))
      used <- list(pointwise = rep(kp, m), simultaneous = rep(ks, m), bonferroni = rep(kb, m))
    } else {
      out <- .np_reg_percentile_e6_materialize(info, sel[seq_len(m)], sel[m + 1L], support[ii, , drop = FALSE],
        matrix(nominal$bands$bonferroni[ii, , , drop = FALSE], m, 2L))
      bands <- out$bands
      z$calibration[[target]] <- list(scores = ss[, m + 1L], empirical.scores = ee[, m + 1L],
        selected = sel[m + 1L], q = sel[m + 1L] / B, marginal.scores = ss[, seq_len(m), drop = FALSE],
        marginal.empirical = ee[, seq_len(m), drop = FALSE], marginal.selected = sel[seq_len(m)],
        marginal.q = sel[seq_len(m)] / B, rank = 38L, J = 39L, attainable = out$attainable,
        support.used = out$support.used, outer.completions = sum(!is.finite(ee) & is.finite(ss)),
        floor = sum(ss[, m + 1L] == 0), support = support[ii, , drop = FALSE], law = 'residual', minimum.content = 0)
      z$ranks[[target]] <- list(kP = out$kp, kS = out$ks, calibration = 'E6 marginal and family rank38; conservative hulls')
      kp <- out$kp; ks <- pmin(kp, out$ks)
      kb <- pmin(nominal$ranks[[target]]$kB, ks)
      if (m == 1L) kp <- kb <- ks
      used <- list(pointwise = kp, simultaneous = ks, bonferroni = kb)
    }
    z$ranks[[target]]$used.indices <- used
    rr <- which(z$resolution$family == target)
    z$resolution[rr, ] <- .np_band_resolution_use_indices(
      z$resolution[rr, , drop = FALSE], used[z$resolution$band[rr]])
    for (b in names(bands)) z$bands[[b]][ii, , ] <- bands[[b]][, , 1L]
  }
  z$center <- p$center; z$se <- p$se; z$expectation <- drop(p$A %*% p$mu)
  z$variance <- rowSums(sweep(p$A, 2L, p$sigma, '*')^2)
  z$map <- p$map; z$budget <- list(bytes = forecast, cap = 256 * 1024^2)
  z$warnings <- character()
  if (any(vapply(z$calibration, function(x) isTRUE(x$support.used), TRUE)))
    z$warnings <- c(z$warnings, 'Full bootstrap-support endpoint used; bounds may be very wide')
  if (any(!vapply(z$calibration, function(x) isTRUE(x$attainable), TRUE)))
    z$warnings <- c(z$warnings, 'Calibration remains unattainable; all bands retained as diagnostic only')
  z
}
