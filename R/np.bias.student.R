# Re-studentization leaves the corrected rows and observed standard errors fixed.
# Pilot/variance providers return at most 16 training rows; no square operator
# or persistent cache is constructed.
.np_bias_student_p0_hat <- function(bws, x, at) {
  .np_plot_with_local_compiled_eval(
    .npreghat_exact_lp_matrix_from_regression_core(bws, txdat=x, exdat=at,
      degree=bws[["degree"]], basis=bws[["basis"]],
      bernstein.basis=bws[["bernstein.basis"]]))
}

.np_bias_student_pilot_owner <- function(pilot) {
  force(pilot)
  function(ii) {
    n <- nrow(pilot$x)
    if (!length(ii) || length(ii)>16L || anyNA(ii) ||
        any(ii<1L | ii>n | ii!=floor(ii)))
      stop("invalid studentization pilot row request", call.=FALSE)
    if (isTRUE(pilot$categorical.only))
      return(.np_bias_categorical_average_rows(pilot$groups,pilot$profile.index[ii],n))
    if (!is.null(pilot$mixed)) return(.np_bias_mixed_p0_owner(pilot$mixed)(ii))
    if (!is.null(pilot$profile)) {
      a <- pilot$profile
      return(.np_bias_profile_pilot_rows(a,a$x[ii,,drop=FALSE],a$h$p0,"pooled P0"))
    }
    if (!is.null(pilot$adaptation)) {
      a <- pilot$adaptation
      return(do.call(rbind, lapply(ii, function(i)
        .np_bias_adaptation_pilot(a, a$x[i,,drop=FALSE], a$geometry$b0,
          a$Q[i,,drop=FALSE])$H)))
    }
    out <- matrix(0,length(ii),n)
    for (group in pilot$groups) {
      rows <- which(ii %in% group)
      if (!length(rows)) next
      H <- .np_bias_student_p0_hat(pilot$bws, pilot$x[group,,drop=FALSE],
        pilot$x[ii[rows],,drop=FALSE])
      ridge <- attr(H,"ridge.used",exact=TRUE)
      if (length(ridge)!=length(rows) || any(!is.finite(ridge)) || any(ridge!=0))
        stop("P0 studentization requires an identified residual pilot", call.=FALSE)
      out[rows,group] <- H
    }
    out
  }
}

.np_bias_student_variance_owner <- function(smoother) {
  force(smoother)
  function(ii) {
    n <- nrow(smoother$x)
    if (!length(ii) || length(ii)>16L || anyNA(ii) ||
        any(ii<1L | ii>n | ii!=floor(ii)))
      stop("invalid studentization variance row request", call.=FALSE)
    out <- matrix(0,length(ii),n)
    for (group in smoother$groups) {
      rows <- which(ii %in% group)
      if (!length(rows)) next
      # Constants cancel in this positive product-kernel normalization.
      # Keep the frozen one-coordinate operation order.
      K <- matrix(1,length(rows),length(group))
      for (j in seq_along(smoother$h))
        K <- K * exp(-.5*(outer(smoother$x[[j]][ii[rows]],
          smoother$x[[j]][group],"-")/smoother$h[j])^2)
      mass <- rowSums(K)
      if (any(!is.finite(K)) || any(!is.finite(mass)) || any(mass<=0))
        stop("invalid studentization variance smoother", call.=FALSE)
      out[rows,group] <- K/mass
    }
    if (!is.null(smoother$q))
      out <- sweep(out,2L,smoother$q,"*")/smoother$denominator[ii]
    out
  }
}

.np_bias_student_prepare <- function(state) {
  if (isTRUE(state$categorical.only)) return(.np_bias_categorical_student_prepare(state))
  if (!is.null(state$mixed)) return(.np_bias_mixed_student_prepare(state))
  if (isTRUE(state$pooled)) return(.np_bias_profile_student_prepare(state))
  pilot <- list(x=state$xcon, groups=state$groups, bws=state$pilot.bws,
    adaptation=if (is.null(state$adaptation)) NULL else
      state$adaptation[c("x","Q","geometry")])
  owner <- .np_bias_student_pilot_owner(pilot)
  # Zero padding preserves the frozen three-response product shape without
  # using a true signal or any other simulation-only input.
  check <- .np_bias_residual_apply(owner,
    cbind(state$y,state$y-state$residual_raw,0))
  response.scale <- max(abs(state$y), .Machine$double.xmin)
  if (max(abs(check$row_variance-state$residual_row_variance))>1e-8 ||
      max(abs(check$residual[,1L]-state$residual)/response.scale)>1e-8)
    stop("studentization residual pilot disagrees with the observed pilot", call.=FALSE)
  smoother <- list(x=state$xcon,groups=state$groups,h=state$variance_h)
  generation <- state$variance
  if (!is.null(state$local.inference)) {
    smoother$q <- state$residual_row_variance
    smoother$denominator <- state$local.inference$denominator
    generation <- state$local.inference$generation
  }
  list(pilot=pilot, q=state$residual_row_variance, rmu=check$residual[,2L],
    root.variance=sqrt(generation),smoother=smoother)
}

.np_bias_student_target <- function(prepared, W, se) {
  A <- .np_bias_right_compose(W^2,
    .np_bias_student_variance_owner(prepared$smoother))
  if (any(!is.finite(A)) || any(A<0))
    stop("invalid studentization variance composition",call.=FALSE)
  active <- rowSums(abs(W))>0
  if (any(!is.finite(se)) || any(se<0) || any(active & se==0) || any(!active & se!=0))
    stop("inconsistent corrected target scale",call.=FALSE)
  list(A=A,active=active)
}

.np_bias_student_slices <- function(size) {
  lapply(seq.int(1L,size,by=256L),function(first)
    seq.int(first,min(size,first+255L)))
}

# The q=0 saturated profile owner uses observed residual wild errors. For
# profile averages their conditional covariance equals the cell-variance map,
# while heterogeneous residual magnitudes retain an estimable replicate scale.
# Existing continuous owners keep their frozen smoothed-variance construction.
.np_bias_student_error_scale <- function(prepared, law) {
  if (identical(law,"mammen") && !is.null(prepared$wild.residual))
    prepared$wild.residual else prepared$root.variance
}

.np_bias_student_errors <- function(root.variance, size, law, innovations=NULL) {
  lapply(.np_bias_student_slices(size),function(ii)
    root.variance * (if (is.null(innovations))
      .np_bias_draws(length(root.variance),length(ii),law) else innovations[,ii,drop=FALSE]))
}

.np_bias_student_chunk <- function(prepared, W, target, size, law, signed=FALSE, innovations=NULL) {
  n <- ncol(W)
  slices <- .np_bias_student_slices(size)
  errors <- .np_bias_student_errors(.np_bias_student_error_scale(prepared,law),size,law,innovations)
  numerator <- lapply(errors,function(ep) W %*% ep)
  residual <- matrix(0,n,size)
  owner <- .np_bias_student_pilot_owner(prepared$pilot)
  for (first in seq.int(1L,n,by=16L)) {
    ii <- seq.int(first,min(n,first+15L))
    H <- owner(ii)
    for (b in seq_along(slices))
      residual[ii,slices[[b]]] <-
        (errors[[b]][ii,,drop=FALSE]-H %*% errors[[b]])/sqrt(prepared$q[ii])
  }
  out <- matrix(0,nrow(W),size)
  for (b in seq_along(slices)) {
    ii <- slices[[b]]
    scale <- sqrt(target$A %*% (residual[,ii,drop=FALSE]+prepared$rmu)^2)
    if (any(!is.finite(scale)) || any(scale<=0 & target$active))
      stop("nonpositive replicate studentization scale",call.=FALSE)
    # Only deterministic zero rows get a harmless divisor: their numerator and
    # returned pivot are exactly zero. Active scales are never floored.
    if (any(!target$active)) scale[!target$active,] <- 1
    out[,ii] <- if(signed) numerator[[b]]/scale else abs(numerator[[b]]/scale)
  }
  out
}

.np_bias_student_draws <- function(prepared, W, se, B, law, signed=FALSE, innovations=NULL) {
  target <- .np_bias_student_target(prepared,W,se)
  .npRmpi_bootstrap_fanout_enabled(n=ncol(W),B=B,chunk.size=1024L,
    what="studentized corrected process")
  tasks <- .npRmpi_bootstrap_chunk_tasks(B,1024L,with.seeds=FALSE)
  if (is.null(innovations)) {
    tasks <- .npRmpi_bootstrap_attach_rng_stream(tasks,function(task)
      invisible(.np_bias_student_errors(.np_bias_student_error_scale(prepared,law),task$bsz,law)),
      normal.draws=identical(law,"gaussian"))
  } else {
    tasks <- lapply(tasks,function(task) {
      ii <- seq.int(task$start,length.out=task$bsz)
      task$innovations <- innovations$values[,ii,drop=FALSE]
      task
    })
    attr(tasks,"rng_final_state") <- innovations$end
  }
  worker <- function(task,prepared,W,target,law,signed) {
    if (is.null(task$innovations)) .npRmpi_bootstrap_use_task_rng(task)
    t(.np_bias_student_chunk(prepared,W,target,task$bsz,law,signed=signed,
      innovations=task$innovations))
  }
  t(.npRmpi_bootstrap_run_fanout(tasks,worker,ncol.out=nrow(W),
    what="studentized corrected process",profile.where="mpi.applyLB:bias-student",
    progress.label="Calibrating studentized process",progress.unit="draws",
    prepared=prepared,W=W,target=target,law=law,signed=signed))
}
