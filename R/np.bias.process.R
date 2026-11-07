# Calibration is estimator-neutral. A family adapter supplies complete rows,
# centers and scales. One observation draw is shared across the whole family.
.np_bias_draws <- function(n, size, law) {
  switch(law,
    gaussian=matrix(rnorm(n*size),n),
    mammen={p<-(sqrt(5)+1)/(2*sqrt(5));
      matrix(ifelse(runif(n*size)<p,(1-sqrt(5))/2,(1+sqrt(5))/2),n)},
    rademacher=matrix(ifelse(runif(n*size)<.5,-1,1),n),
    stop("unknown bias multiplier law",call.=FALSE))
}

# The same bounded innovation/chunk contract is used by the MPI transport.
.np_bias_process_innovations <- function(n, first, last, B, law, blocklen) {
  if (law == "inid") return(rmultinom(last-first+1L,n,rep(1/n,n))-1)
  if (law %in% c("fixed","geom"))
    return(.np_block_counts_drawer(n,B,blocklen,law)(first,last)-1)
  .np_bias_draws(n,last-first+1L,law)
}

.np_bias_process_chunk <- function(F, first, last, B, law, blocklen) {
  F %*% .np_bias_process_innovations(ncol(F),first,last,B,law,blocklen)
}

# Box-Muller has a cached deviate outside .Random.seed. Retain actual values
# for bounded Gaussian regression jobs; no seed reset can replay that cache.
# Other laws/kinds retain their existing path. Refuse uncached Box-Muller
# here using resolved dimensions, before internal draws or seed-only replay.
.np_bias_gaussian_cache <- function(n, B, law, kind=RNGkind()) {
  builtin <- c("Wichmann-Hill","Marsaglia-Multicarry","Super-Duper",
    "Mersenne-Twister","Knuth-TAOCP","Knuth-TAOCP-2002","L'Ecuyer-CMRG")
  if (!identical(law,"gaussian") || !kind[1L] %in% builtin ||
      !identical(kind[2L],"Box-Muller")) return(NULL)
  dims <- c(n,B)
  if (!is.numeric(dims) || length(dims)!=2L || anyNA(dims) ||
      any(!is.finite(dims)) || any(dims<1 | dims!=floor(dims)) ||
      any(dims>.Machine$integer.max))
    stop("invalid Gaussian replay dimensions",call.=FALSE)
  if (8*as.double(n)*B>16*1024^2)
    stop("unsafe RNG replay: resolved Gaussian regression exceeds the 16 MiB Box-Muller innovation cache; explicitly select a qualified normal generator such as RNGkind(normal.kind = 'Inversion') before plotting",call.=FALSE)
  values <- .np_bias_draws(n,B,"gaussian")
  list(values=values,end=get(".Random.seed",.GlobalEnv,inherits=FALSE))
}

.np_bias_process_draws <- function(F, B, law, blocklen, innovations=NULL) {
  n <- ncol(F)
  .npRmpi_bootstrap_fanout_enabled(n=n,B=B,chunk.size=32L,what="corrected process")
  tasks <- .npRmpi_bootstrap_chunk_tasks(B,32L,with.seeds=FALSE)
  if (is.null(innovations)) {
    tasks <- .npRmpi_bootstrap_attach_rng_stream(tasks,function(task)
      invisible(.np_bias_process_innovations(n,task$start,task$start+task$bsz-1L,B,law,blocklen)),
      normal.draws=identical(law,"gaussian") || (identical(law,"geom") && blocklen>1L))
  } else {
    tasks <- lapply(tasks,function(task) {
      ii <- seq.int(task$start,length.out=task$bsz)
      task$innovations <- innovations$values[,ii,drop=FALSE]
      task
    })
    attr(tasks,"rng_final_state") <- innovations$end
  }
  worker <- function(task,factors,B,law,blocklen) {
    if (!is.null(task$innovations)) return(t(factors %*% task$innovations))
    .npRmpi_bootstrap_use_task_rng(task)
    t(.np_bias_process_chunk(factors,task$start,task$start+task$bsz-1L,B,law,blocklen))
  }
  t(.npRmpi_bootstrap_run_fanout(tasks,worker,ncol.out=nrow(F),
    what="corrected process",profile.where="mpi.applyLB:bias-process",
    progress.label="Calibrating corrected process",factors=F,B=B,law=law,blocklen=blocklen))
}

.np_bias_block_selection <- function(residual, method, blocklen=NULL) {
  column <- if(method=="geom") "BstarSB" else "BstarCB"
  if(!is.null(blocklen)) {
    if(!is.numeric(blocklen)||length(blocklen)!=1L||!is.finite(blocklen)||blocklen<1||
       blocklen>length(residual)||blocklen!=floor(blocklen))stop("invalid explicit block length",call.=FALSE)
    return(list(length=as.integer(blocklen),source="override",column=column,method=method))
  }
  if(any(!is.finite(residual))||var(residual)<=0)stop("block selection requires nonconstant finite residuals",call.=FALSE)
  selected <- b.star(matrix(residual,ncol=1L),round=TRUE)
  value <- unname(selected[1L,column])
  if(!is.finite(value)||value<1||value>length(residual)||value!=floor(value))
    stop("b.star returned an invalid block length",call.=FALSE)
  list(length=as.integer(value),source="raw P0 residuals",column=column,
    method=method,selector="b.star",both=selected)
}

# Retain exact marginal order statistics and one draw maximum per family.
# The full grid-by-draw matrix is never retained by inference.
.np_bias_limits_prepare <- function(groups, alpha, gaussian, B) {
  index <- split(seq_along(groups),groups)
  list(index=index,alpha=alpha,gaussian=gaussian,
    maximum=matrix(0,length(index),B),
    point=matrix(0,length(groups),length(alpha)),
    bonf=matrix(0,length(groups),length(alpha)))
}

.np_bias_limits_update <- function(state, absolute, rows) {
  if(nrow(absolute)!=length(rows) || ncol(absolute)!=ncol(state$maximum) ||
     any(!is.finite(absolute)) || any(absolute<0))
    stop("invalid corrected process draws",call.=FALSE)
  for(g in seq_along(state$index)) {
    local <- which(rows %in% state$index[[g]])
    if(!length(local))next
    maximum <- state$maximum[g,]
    for(j in local)maximum <- pmax(maximum,absolute[j,])
    state$maximum[g,] <- maximum
    if(!state$gaussian)for(j in local) {
      q <- unname(quantile(absolute[j,],c(1-state$alpha,1-state$alpha/length(state$index[[g]])),type=1L))
      state$point[rows[j],] <- q[seq_along(state$alpha)]
      state$bonf[rows[j],] <- q[length(state$alpha)+seq_along(state$alpha)]
    }
  }
  state
}

.np_bias_limits_finish <- function(state) {
  alpha <- state$alpha;index <- state$index
  limits <- setNames(lapply(c("pointwise","simultaneous","bonferroni"),function(b)
    array(0,c(nrow(state$point),2L,length(alpha)),dimnames=list(NULL,c("lower","upper"),as.character(alpha)))),
    c("pointwise","simultaneous","bonferroni"))
  raw.max <- matrix(0,length(alpha),length(index),dimnames=list(as.character(alpha),names(index)))
  for(g in seq_along(index)) {
    ii <- index[[g]];m <- length(ii)
    raw.max[,g] <- unname(quantile(state$maximum[g,],1-alpha,type=1L))
    if(state$gaussian) {
      point <- matrix(rep(qnorm(alpha/2,lower.tail=FALSE),each=m),m,length(alpha))
      bonf <- matrix(rep(qnorm(alpha/(2*m),lower.tail=FALSE),each=m),m,length(alpha))
      simultaneous <- pmin(bonf,pmax(point,matrix(rep(raw.max[,g],each=m),m,length(alpha))))
    } else {
      point <- state$point[ii,,drop=FALSE];bonf <- state$bonf[ii,,drop=FALSE]
      simultaneous <- matrix(rep(raw.max[,g],each=m),m,length(alpha))
      if(any(point>simultaneous))stop("maximum marginal-order invariant failed",call.=FALSE)
      bonf <- pmax(bonf,simultaneous)
    }
    for(band in names(limits)) {
      critical <- switch(band,pointwise=point,simultaneous=simultaneous,bonferroni=bonf)
      limits[[band]][ii,1L,] <- -critical
      limits[[band]][ii,2L,] <- critical
    }
  }
  list(limits=limits,raw_max_critical=raw.max)
}

.np_bias_limits <- function(absolute,groups,alpha,gaussian=FALSE) {
  state <- .np_bias_limits_prepare(groups,alpha,gaussian,ncol(absolute))
  .np_bias_limits_finish(.np_bias_limits_update(state,absolute,seq_len(nrow(absolute))))
}

# Replaying a draw stream for bounded target batches does not consume a new
# independent stream. Advance the caller once, by one complete draw sequence.
# Repeated target/family draws require all normal state, not only .Random.seed.
# R encodes the normal generator in the seed header's hundreds field (2 is
# Box-Muller). Inspect the saved stream: it may differ from the caller's kind.
# User callbacks may keep hidden state even when their seed payload is nonempty.
.np_bias_custom_replay_check <- function(custom.uniform, custom.normal, normal) {
  if (isTRUE(custom.uniform) || (isTRUE(normal) && isTRUE(custom.normal)))
    stop("unsafe RNG replay: user-supplied generator state cannot be certified by a seed vector; explicitly select a qualified built-in generator before this replay", call. = FALSE)
  invisible(NULL)
}

.np_bias_normal_replay_check <- function(law, blocklen, replay, retained = FALSE,
                                          seed = NULL) {
  if (!isTRUE(replay) || isTRUE(retained)) return(invisible(NULL))
  normal <- identical(law, "gaussian") ||
    (identical(law, "geom") && isTRUE(blocklen > 1))
  if (!normal) return(invisible(NULL))
  box.muller <- if (is.null(seed)) identical(RNGkind()[2L], "Box-Muller") else
    isTRUE((seed[1L] %/% 100L) %% 100L == 2L)
  if (box.muller)
    stop("unsafe RNG replay: seed-only replay cannot reproduce Box-Muller normal state; explicitly select a qualified normal generator such as RNGkind(normal.kind = 'Inversion') before plotting", call. = FALSE)
  invisible(NULL)
}

# Known missing state and uncertifiable custom state have distinct refusals.
# A nonempty custom payload is not a certificate of RNG replayability.
.np_bias_replay_check <- function(law, blocklen, replay, retained = FALSE,
                                seed = NULL) {
  if (!isTRUE(replay) || isTRUE(retained)) return(invisible(NULL))
  .npRmpi_replay_context_preflight()
  if (is.null(seed)) {
    kind <- RNGkind()
    absent <- .np_plot_replay_uniform_state_absent(kind)
    custom.normal <- identical(kind[2L], "user-supplied")
    custom.uniform <- identical(kind[1L], "user-supplied")
  } else {
    absent <- is.integer(seed) && length(seed) == 1L &&
      isTRUE(seed[1L] %% 100L == 5L)
    custom.normal <- isTRUE((seed[1L] %/% 100L) %% 100L == 3L)
    custom.uniform <- isTRUE(seed[1L] %% 100L == 5L)
  }
  # User normal callbacks need not consume the selected uniform generator.
  # Their hidden-state contract is separate from this known-missing-state test.
  if (absent && !(identical(law, "gaussian") && custom.normal))
    stop("unsafe RNG replay: repeated draws require custom uniform state absent from .Random.seed; explicitly select a qualified built-in generator before plotting", call. = FALSE)
  .np_bias_normal_replay_check(law, blocklen, replay, retained, seed)
  .np_bias_custom_replay_check(custom.uniform, custom.normal,
    identical(law, "gaussian") || (identical(law, "geom") && isTRUE(blocklen > 1)))
}

.np_bias_stream <- function(state, code) {
  previous <- get(".Random.seed",.GlobalEnv,inherits=FALSE)
  on.exit(assign(".Random.seed",previous,.GlobalEnv))
  assign(".Random.seed",state,.GlobalEnv)
  value <- force(code)
  list(value=value,end=get(".Random.seed",.GlobalEnv,inherits=FALSE))
}

# Partition both axes without changing the target-major family order or draws.
.np_bias_target_blocks <- function(n.grid, n.targets, anchors = 25L, rows = 500L) {
  if (any(!is.finite(c(n.grid,n.targets,anchors,rows))) ||
      any(c(n.grid,n.targets,anchors,rows) < 1) ||
      any(c(n.grid,n.targets,anchors,rows) != floor(c(n.grid,n.targets,anchors,rows))) ||
      n.grid > .Machine$integer.max / n.targets)
    stop("invalid or unrepresentable bias target grid", call.=FALSE)
  count <- min(n.targets, rows)
  grid.block <- min(anchors, rows %/% count)
  out <- list()
  for (first in seq.int(1L,n.grid,by=grid.block)) {
    jj <- seq.int(first,min(n.grid,first+grid.block-1L))
    for (target in seq.int(1L,n.targets,by=count)) {
      tt <- seq.int(target,min(n.targets,target+count-1L))
      out[[length(out)+1L]] <- list(grid=jj,
        rows=as.vector(outer(jj,(as.double(tt)-1)*n.grid,"+")))
    }
  }
  out
}

.np_bias_inference <- function(state, grid, targets="mean", alpha=.05,
                               B=9999L, law="gaussian", blocklen=NULL, fit=NULL,
                               studentize=TRUE, main.capture=NULL,
                               resolution.bands=c("pointwise","simultaneous","bonferroni"),
                               resolution.warn=TRUE, calibration="absolute", refinement="none") {
  if (!is.character(refinement) || length(refinement)!=1L || is.na(refinement) ||
      !refinement %in% c("none","pilot-block"))
    stop("unknown regression provider refinement",call.=FALSE)
  pilot.block <- identical(refinement,"pilot-block")
  if (pilot.block && !law %in% c("fixed","geom"))
    stop("pilot-block refinement requires fixed or geom",call.=FALSE)
  tail.calibration <- calibration %in% c("tail-envelope","separate-tails")
  if(tail.calibration && !law %in% c("mammen","rademacher","inid","fixed","geom"))
    stop("private separate-tail regression requires a supported bootstrap law",call.=FALSE)
  rank.calibration <- identical(calibration,"rank") || tail.calibration
  if(!rank.calibration && !identical(calibration,"absolute"))
    stop("unknown private regression calibration",call.=FALSE)
  if(identical(calibration,"rank") && (!identical(law,"mammen") || !isTRUE(studentize)))
    stop("private signed-rank regression requires studentized Mammen",call.=FALSE)
  if ((!is.null(state$local.inference) || !is.null(state$main.capture) || !is.null(main.capture)) && !law %in% c("gaussian","mammen"))
    stop("captured local calibration requires Gaussian or Mammen; prepare the state for the selected law",call.=FALSE)
  if ((!is.null(state$local.inference) || !is.null(state$main.capture)) && is.null(main.capture))
    main.capture <- .np_bias_capture_prepare(state,fit,grid,targets)
  if (!is.null(main.capture)) .np_bias_capture_validate(main.capture,state,grid,targets)
  if (isTRUE(state$categorical.only) && !law %in% c("gaussian","mammen","inid","rademacher","fixed","geom"))
    stop("categorical-only bias correction does not support this calibration law",call.=FALSE)
  if (!is.null(state$mixed) && !law %in% c("gaussian","mammen"))
    stop("mixed local-to-global adaptation supports Gaussian or Mammen calibration; use adaptation='none' for incumbent laws",call.=FALSE)
  higher <- !is.null(state$gradient.order) && any(vapply(targets,function(t)
    startsWith(t,"d") && .np_bias_derivative_order(state,as.integer(substring(t,2L)))>1L,logical(1)))
  if (higher && !law %in% c("gaussian","mammen","inid","rademacher","fixed","geom"))
    stop("higher-order bias correction does not support this calibration law",call.=FALSE)
  if (isTRUE(state$pooled) &&
      !law %in% c("gaussian","mammen","inid","rademacher","fixed","geom"))
    stop("pooled mixed-data bias correction currently supports Gaussian, Mammen, Rademacher, fixed, geom or inid calibration only",call.=FALSE)
  if(!is.numeric(B)||length(B)!=1L||!is.finite(B)||B<1||B!=floor(B)||B>.Machine$integer.max)
    stop("bias calibration draws must be a positive integer",call.=FALSE)
  if(B<2 && law %in% c("inid","fixed","geom") &&
     !(isTRUE(state$pooled) && law %in% c("fixed","geom")))
    stop("score-bootstrap variance requires at least two draws",call.=FALSE)
  if(!length(alpha)||any(!is.finite(alpha))||any(alpha<=0|alpha>=.5)||anyDuplicated(alpha))
    stop("alpha must contain distinct values in (0, 0.5)",call.=FALSE)
  if(!law%in%c("gaussian","mammen","rademacher","inid","fixed","geom"))stop("unsupported corrected inference law",call.=FALSE)
  if(!identical(names(grid),names(state$x))||!nrow(grid)||anyNA(grid))stop("invalid bias evaluation grid",call.=FALSE)
  if(!is.character(targets)||!length(targets)||anyNA(targets)||anyDuplicated(targets))
    stop("invalid or duplicate bias target family",call.=FALSE)
  if(rank.calibration &&
     8*as.double(B)*(length(targets)+10*min(8L,length(targets)))>128*1024^2)
    stop("rank workspace exceeds declared byte budget before allocation",call.=FALSE)
  blocks <- .np_bias_target_blocks(nrow(grid),length(targets))
  map <- expand.grid(row=seq_len(nrow(grid)),target=targets,stringsAsFactors=FALSE)
  M <- nrow(map);B <- as.integer(B)
  center <- original <- se <- score.se <- numeric(M)
  student <- isTRUE(studentize) && (law %in% c("gaussian","mammen") ||
    (isTRUE(state$pooled) && law=="rademacher"))
  score.calibration <- isTRUE(state$pooled) && law %in% c("fixed","geom")
  calibration.state <- .np_bias_limits_prepare(map$target,alpha,law=="gaussian" && !student,B)
  prepared <- NULL
  regularization <- adaptation <- list()
  block.refinement <- if(pilot.block)list() else NULL
  adapted.rows <- if (!is.null(state$mixed)) .np_bias_mixed_capture(state,fit,grid,targets) else if (!is.null(state$adaptation))
    .np_bias_adaptation_capture(state,fit,grid,targets) else NULL
  score <- law%in%c("inid","fixed","geom")
  selection <- if(law%in%c("fixed","geom")) .np_bias_block_selection(state$residual_raw,law,blocklen) else NULL
  if(!exists(".Random.seed",.GlobalEnv,inherits=FALSE))runif(1L)
  rng <- get(".Random.seed",.GlobalEnv,inherits=FALSE);rng.kind <- RNGkind();end <- rng
  innovations <- .np_bias_gaussian_cache(nrow(state$x),B,law,rng.kind)
  rank.results <- if(rank.calibration)vector("list",length(targets)) else NULL
  families <- if(rank.calibration)split(seq_along(targets),ceiling(seq_along(targets)/8L)) else list(1L)
  for(family in families) {
   if(rank.calibration) {
    rank.states <- lapply(family,function(g)
      (if(tail.calibration).np_band_tail_prepare else .np_band_rank_prepare)(B,map[map$target==targets[g],,drop=FALSE],alpha,
        byte.limit=128*1024^2/length(family)))
    family.rows <- which(map$target %in% targets[family])
    # One capped signed tile can span several existing 25-anchor row blocks.
    # If it fits, retain this tile only and avoid recomputing the process.
    # The cap is independent of grid size; larger families use exact replay.
    retain.tile <- !tail.calibration && length(family.rows)<=128L &&
      16*as.double(length(family.rows))*B+80*as.double(B)*length(family)<=128*1024^2
    signed.tile <- if(retain.tile)matrix(0,length(family.rows),B) else NULL
    tile.rows <- floor((128*1024^2-(if(tail.calibration)96 else 72)*as.double(B)*length(family)-
      if(tail.calibration)48*M*length(alpha) else 0)/(8*B))
    if(tile.rows<length(family))stop("rank workspace exceeds declared byte budget before allocation",call.=FALSE)
    family.blocks <- .np_bias_target_blocks(nrow(grid),length(family),
      anchors=min(25L,tile.rows %/% length(family)),rows=min(500L,tile.rows))
    family.blocks <- lapply(family.blocks,function(b) {b$rows<-family.rows[b$rows];b})
   } else {family.blocks <- blocks;retain.tile <- FALSE}
   .np_bias_replay_check(law, if (is.null(selection)) NULL else selection$length,
     length(family.blocks) > 1L || length(families) > 1L ||
       (rank.calibration && !tail.calibration && !retain.tile),
     retained = !is.null(innovations), seed = rng)
   for(pass in seq_len(if(rank.calibration && !tail.calibration && !retain.tile)2L else 1L)) {
    if(rank.calibration && pass==2L)rank.states <- lapply(rank.states,.np_band_rank_cutoff)
    for(batch in family.blocks) {
    jj <- batch$grid;rr <- batch$rows
    local.map <- map[rr,,drop=FALSE];local.map$row <- match(local.map$row,jj)
    value <- if (!is.null(main.capture)) main.capture$rows(jj,local.map) else if (is.null(adapted.rows)) .np_bias_regression_rows(state,grid[jj,,drop=FALSE],local.map) else
      adapted.rows(jj,local.map)
    if (pass==1L && !is.null(value$adaptation)) adaptation[[length(adaptation)+1L]] <-
      list(rows=rr,grid=jj,details=value$adaptation)
    diagnostic <- value$main.regularization
    if(law!="gaussian" && any(diagnostic$ridge>0) && is.null(state$local.inference) && is.null(main.capture))
      stop("bootstrap calibration for regularized generalized-NN main rows remains unqualified; Gaussian inference is available with the qualifications documented in np_bias_control",call.=FALSE)
    diagnostic$row <- jj[diagnostic$row]
    if(pass==1L)regularization[[length(regularization)+1L]] <- diagnostic
    center[rr] <- value$center;original[rr] <- value$original_center;se[rr] <- value$se
    divisor <- ifelse(value$se==0,1,value$se)
    if (student) {
      if (is.null(prepared)) prepared <- .np_bias_student_prepare(state)
      part <- if(rank.calibration) .np_bias_stream(rng,
        .np_bias_student_draws(prepared,value$W,value$se,B,law,signed=TRUE,innovations=innovations)) else
        .np_bias_stream(rng,.np_bias_student_draws(prepared,value$W,value$se,B,law,innovations=innovations))
    } else if (score.calibration) {
      F <- sweep(value$W,2L,state$residual,"*")
      F <- F-rowMeans(F)
      part <- if(rank.calibration) .np_bias_stream(rng,
        .np_bias_score_student_draws(F,B,law,selection$length,value$se>0,studentize,signed=TRUE)) else
        .np_bias_stream(rng,
          .np_bias_score_student_draws(F,B,law,selection$length,value$se>0,studentize))
      se[rr] <- score.se[rr] <- part$value$se
      part$value <- if(rank.calibration) part$value$signed else part$value$absolute
    } else {
      F <- if(score)sweep(value$W,2L,state$residual,"*") else value$factor
      if(score)F <- F-rowMeans(F)
      part <- .np_bias_stream(rng,
        .np_bias_process_draws(F,B,law,if(is.null(selection))NULL else selection$length,innovations=innovations))
    }
    .np_plot_progress_notify(stage="Constructing bias-corrected intervals",
      done=max(jj),total=nrow(grid),unit="grid")
    end <- part$end
    if(score && !score.calibration) {
      s <- apply(part$value,1L,sd)
      if(any(!is.finite(s))||any(s<=0 & value$se>0))stop("nonpositive score-bootstrap scale",call.=FALSE)
      score.se[rr] <- s
      divisor <- ifelse(s==0,1,s)
      if(law!="inid")se[rr] <- s
      if (pilot.block) {
        anchored <- .np_bias_pilot_block_scale(F,value$se,selection$length,law)
        se[rr] <- anchored$se
        if (pass==1L) block.refinement[[length(block.refinement)+1L]] <- c(list(rows=rr),anchored)
      }
    }
    if(rank.calibration) {
      if(!student && !score.calibration) part$value <- part$value/divisor
      # A zero target may share nonzero family cutoffs, but never nonzero draws.
      for(zero in which(se[rr]==0)) if(any(part$value[zero,]!=0))
        stop("zero target scale with nonzero signed process draws",call.=FALSE)
      if(retain.tile)signed.tile[match(rr,family.rows),] <- part$value
      for(g in seq_along(family)) {
        local <- which(map$target[rr]==targets[family[g]])
        if(length(local)) rank.states[[g]] <- if(tail.calibration)
          .np_band_tail_update(rank.states[[g]],part$value[local,,drop=FALSE],
            map[rr[local],,drop=FALSE],seq_len(B)) else
          .np_band_rank_update(rank.states[[g]],part$value[local,,drop=FALSE],
            map[rr[local],,drop=FALSE],seq_len(B),replay=pass==2L)
      }
      if(pass==1L)calibration.state <- .np_bias_limits_update(calibration.state,abs(part$value),rr)
    } else calibration.state <- .np_bias_limits_update(calibration.state,
      if(student || score.calibration)part$value else abs(part$value/divisor),rr)
    }
   }
   if(rank.calibration && retain.tile) {
    rank.states <- lapply(rank.states,.np_band_rank_cutoff)
    for(g in seq_along(family)) {
      local <- which(map$target[family.rows]==targets[family[g]])
      rank.states[[g]] <- .np_band_rank_update(rank.states[[g]],
        signed.tile[local,,drop=FALSE],map[family.rows[local],,drop=FALSE],seq_len(B),replay=TRUE)
    }
    signed.tile <- NULL
   }
   if(rank.calibration)for(g in seq_along(family))rank.results[[family[g]]] <-
      (if(tail.calibration).np_band_tail_finish else .np_band_rank_finish)(rank.states[[g]])
  }
  assign(".Random.seed",end,.GlobalEnv)
  calibrated <- .np_bias_limits_finish(calibration.state)
  absolute.reference <- NULL
  if(rank.calibration) {
    absolute.reference <- calibrated
    root.limits <- calibrated$limits
    for(g in seq_along(targets))for(band in names(root.limits))
      root.limits[[band]][map$target==targets[g],,] <- rank.results[[g]]$limits[[band]]
    calibrated <- if(identical(calibration,"tail-envelope"))
      c(.np_band_tail_protect(root.limits,absolute.reference$limits),list(families=rank.results)) else
      list(limits=root.limits,families=rank.results)
  }
  bands <- if(rank.calibration).np_band_rank_invert(center,se,calibrated$limits,deterministic=se==0) else lapply(calibrated$limits,function(x) {
    for(k in seq_along(alpha))for(side in 1:2)x[,side,k] <- center+se*x[,side,k]
    x
  })
  resolution <- if(length(resolution.bands))
    .np_band_resolution(B,map$target,alpha,se>0,resolution.bands,
      signed=rank.calibration,simultaneous.index=if(rank.calibration) {
        k <- do.call(rbind,lapply(rank.results,function(x)x$kS));rownames(k)<-targets
        k[names(split(seq_along(map$target),map$target)),,drop=FALSE]
      } else NULL) else NULL
  if(!is.null(resolution) && law=="gaussian" && !student) {
    analytic <- resolution$band %in% c("pointwise","bonferroni")
    resolution$nominal_tail_draws[analytic] <- NA_real_
    resolution$order_index[analytic] <- NA_real_
    resolution$resolution_limited[analytic] <- FALSE
    resolution$calibration[analytic] <- "analytic Gaussian"
  }
  if(isTRUE(resolution.warn)) .np_band_resolution_warn(resolution)
  result <- list(resolution=resolution,construction="iterated",version=if(score.calibration)3L else if(student)2L else 1L,estimate=original,center=center,
    correction=original-center,se=se,bands=bands,map=map,grid=grid,
    settings=state$settings,law=law,block.selection=selection,draws=B,alpha=alpha,
    main.regularization=do.call(rbind,regularization),
    adaptation=if(length(adaptation))adaptation else NULL,
    covariance=list(operator="complete corrected observation rows",
      variance=if(score.calibration)"exact conditional selected-block score variance" else
        if(law%in%c("fixed","geom"))"empirical residual-score process variance" else
        if(isTRUE(state$categorical.only))"within-profile residual variance" else
        if(!is.null(state$local.inference))"pooled normalized P0 residual variance" else "smoothed P0 residual variance",
      studentization=if(score.calibration) {
        if(studentize)"replicate-specific selected-block score SD" else "observed selected-block score SD"
      } else if(isTRUE(state$categorical.only) && student)"replicate-specific within-profile residual SD" else
        if(score)"empirical residual-score process SD" else
        if(isTRUE(state$categorical.only))"within-profile residual SD" else
        if(student)"replicate-specific smoothed P0 residual SD" else "smoothed complete-row SD",
      calibration=if(score.calibration) {
        if(studentize)"conditional score bootstrap-t process" else "conditional score process"
      } else if(student)"residual-pilot studentized multiplier process" else
        if(score)"residual score process" else "fixed-factor multiplier process"),
    studentization=if(score.calibration && studentize)list(construction="conditional-block-score-v1",
      draw.block=32L) else if(student)list(construction=if(isTRUE(state$categorical.only))
      "finite-profile-studentized-v1" else if(!is.null(state$local.inference))
      "captured-pooled-positive-v1" else "residual-pilot-v1",
      draw.block=256L,reuse.block=1024L,training.block=16L) else NULL,
    calibration=calibrated,score.se=if(score)score.se else NULL,
    rng=list(kind=rng.kind,start=rng,end=end,convention="common observation draws across all target blocks"),
    scope="finite grid within each scalar-target family",owner="bounded corrected observation factors")
  if (pilot.block) {
    result$provider.refinement <- list(method=refinement,
      owner=if(score.calibration)"unchanged pooled score studentizer" else "pilot times selected-block to iid covariance ratio",
      blocks=block.refinement)
    if (!score.calibration) result$covariance$variance <- "pilot-anchored exact selected-block relative covariance"
  }
  if(rank.calibration) {
    result$absolute.reference <- absolute.reference
    result$calibration.method <- calibration
    if(any(se==0)) result$deterministic.targets <- which(se==0)
  }
  result
}

.np_bias_plot_errors <- function(state, grid, target, method, bootstrap, wild,
                                 alpha, B, gaussian.draws, blocklen, band, fit=NULL,
                                 resolution.context=NULL, calibration="default") {
  law <- if(method=="asymptotic")"gaussian" else if(bootstrap=="wild")wild else bootstrap
  effective <- if(calibration=="default")"absolute" else
    .np_plot_calibration_method(calibration, "regression")
  refinement <- if(law %in% c("fixed","geom") && band!="pmzsd")"pilot-block" else "none"
  inference <- .np_plot_corrected_activity_run("Constructing bias-corrected intervals",
    if (identical(effective, "smoothed-percentile"))
      .np_reg_percentile_plot(state, fit, grid, target, B,
        if (band=="all") c("pointwise","simultaneous","bonferroni") else band,
        resolution.context=resolution.context) else .np_bias_inference(state,grid,target,alpha,
      if(method=="asymptotic")gaussian.draws else B,law,blocklen,fit,
      studentize=band!="pmzsd",
      resolution.bands=if(band=="pmzsd")character() else if(band=="all")
        c("pointwise","simultaneous","bonferroni") else band,
      resolution.warn=is.null(resolution.context),calibration=effective,refinement=refinement))
  inference$calibration.requested <- calibration
  inference$calibration.method <- effective
  if (effective %in% c("separate-tails","tail-envelope")) {
    inference$root.orientation <- "estimate minus truth; root-inverted endpoints"
    inference$resolution$calibration <- effective
  }
  if(!is.null(resolution.context) && !is.null(inference$resolution))
    resolution.context$rows[[length(resolution.context$rows)+1L]] <- inference$resolution
  all <- lapply(inference$bands,function(x)cbind(inference$center-x[,1,1],x[,2,1]-inference$center))
  width <- if(band=="pmzsd")cbind(inference$se,inference$se)*qnorm(alpha/2,lower.tail=FALSE) else
    all[[if(band=="all")"pointwise" else band]]
  if(is.null(width))stop("unsupported bias interval type",call.=FALSE)
  list(boot.err=cbind(width,inference$center),boot.all.err=if(band=="all")all else NULL,
    inference=inference,bxp=list())
}
