# This function implements an individual test of significance for both
# discrete (Racine, hart, Li, 2006, ER) and continuous variables
# (Racine, 1997, JBES). It accepts a data frame for explanatory data
# (mixed datatypes allowed), a vector for y for a regression model, an
# npregbw object, and a set of indices for the columns of X for which
# the test is to be run (default = all).

if (getRversion() >= "2.15.1")
  utils::globalVariables(".npsig_worker")

.npRmpi_npsig_extract_xy_from_npreg <- function(obj) {
  if (is.null(obj$bws$formula) || is.null(obj$bws$call))
    stop("unable to extract xdat/ydat from npreg object")

  tt <- terms(obj$bws$formula)
  m <- match(c("formula", "data", "subset", "na.action"),
             names(obj$bws$call), nomatch = 0)
  tmf <- obj$bws$call[c(1, m)]
  tmf[[1]] <- as.name("model.frame")
  tmf[["formula"]] <- tt
  mf.args <- as.list(tmf)[-1L]
  tmf <- do.call(stats::model.frame, mf.args, envir = environment(tt))

  ydat <- model.response(tmf)
  xdat <- tmf[, attr(attr(tmf, "terms"), "term.labels"), drop = FALSE]
  list(xdat = xdat, ydat = ydat)
}

.npRmpi_npsig_extract_xy_from_bws <- function(obj) {
  if (!is.null(obj$call)) {
    call.names <- names(obj$call)
    if (!is.null(call.names) &&
        any(call.names == "xdat") &&
        any(call.names == "ydat")) {
      xdat <- .np_eval_bws_call_arg(obj, "xdat")
      ydat <- .np_eval_bws_call_arg(obj, "ydat")
      return(list(xdat = xdat, ydat = ydat))
    }
  }

  if (!is.null(obj$formula) && !is.null(obj$call)) {
    tt <- terms(obj$formula)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(obj$call), nomatch = 0)
    tmf <- obj$call[c(1, m)]
    tmf[[1]] <- as.name("model.frame")
    tmf[["formula"]] <- tt
    mf.args <- as.list(tmf)[-1L]
    tmf <- do.call(stats::model.frame, mf.args, envir = environment(tt))

    ydat <- model.response(tmf)
    xdat <- tmf[, attr(attr(tmf, "terms"), "term.labels"), drop = FALSE]
    return(list(xdat = xdat, ydat = ydat))
  }

  stop("unable to extract xdat/ydat from bandwidth object")
}

.npRmpi_npsig_validate_index <- function(index, xdat) {
  if(anyNA(index)) stop("index must not contain missing values")
  if(any(index < 1 | index > NCOL(xdat), na.rm = TRUE)) stop(paste("invalid index provided: index entries must lie between 1 and ",NCOL(xdat),sep=""))
  if(length(unique(index)) < length(index)) stop("index contains repeated values (must be unique)")
  invisible(TRUE)
}

.npRmpi_with_local_regression <- function(expr) {
  old.disable <- getOption("npRmpi.autodispatch.disable", FALSE)
  old.ctx <- getOption("npRmpi.autodispatch.context", FALSE)
  old.local <- getOption("npRmpi.local.regression.mode", FALSE)
  options(npRmpi.autodispatch.disable = TRUE)
  options(npRmpi.autodispatch.context = TRUE)
  options(npRmpi.local.regression.mode = TRUE)
  on.exit(options(npRmpi.autodispatch.disable = old.disable), add = TRUE)
  on.exit(options(npRmpi.autodispatch.context = old.ctx), add = TRUE)
  on.exit(options(npRmpi.local.regression.mode = old.local), add = TRUE)
  old.mode <- .Call("C_np_set_local_regression_mode", TRUE, PACKAGE = "npRmpi")
  on.exit(.Call("C_np_set_local_regression_mode", old.mode, PACKAGE = "npRmpi"), add = TRUE)
  force(expr)
}

.npRmpi_npsig_npreg_local <- function(...) {
  .npRmpi_with_local_regression(npreg(...))
}

.npRmpi_npsig_do_local <- function(extra.args = NULL, ...) {
  args <- c(list(...), if (length(extra.args)) extra.args else NULL)
  do.call(.npRmpi_npsig_npreg_local, args)
}

.npRmpi_npsig_do_leaf <- function(fun, extra.args = NULL, ...) {
  args <- c(list(...), if (length(extra.args)) extra.args else NULL)
  .npRmpi_autodispatch_untag(do.call(fun, args))
}

.npRmpi_npsig_npreg_leaf <- function(extra.args = NULL, ...) {
  .npRmpi_npsig_do_leaf(npreg, extra.args = extra.args, ...)
}

.npRmpi_npsig_collective_context <- function() {
  isTRUE(.npRmpi_autodispatch_called_from_bcast())
}

.npRmpi_npsig_bootstrap_seed_plan <- function(num.obs,
                                              boot.num,
                                              boot.method,
                                              draw.wild.mult,
                                              a,
                                              b,
                                              p.a) {
  seeds <- vector("list", boot.num)
  for (i.star in seq_len(boot.num)) {
    seeds[[i.star]] <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (boot.method == "iid" || boot.method == "pairwise") {
      sample.int(num.obs, replace = TRUE)
    } else if (boot.method == "wild") {
      draw.wild.mult(num.obs, a, b, p.a)
    } else if (boot.method == "wild-rademacher") {
      draw.wild.mult(num.obs, -1, 1, p.a)
    } else {
      stop(sprintf("unsupported bootstrap method '%s'", boot.method), call. = FALSE)
    }
  }
  seeds
}

.npRmpi_npsig_gather_rank_chunks <- function(gathered, size) {
  size <- as.integer(size)[1L]
  if (is.na(size) || size < 1L)
    stop("invalid MPI gather size")

  if (is.matrix(gathered)) {
    if (!identical(ncol(gathered), size))
      stop("npsigtest MPI gather returned malformed matrix output", call. = FALSE)
    return(lapply(seq_len(size), function(j) gathered[, j]))
  }

  if (is.array(gathered)) {
    dims <- dim(gathered)
    if (length(dims) < 2L || !identical(dims[[length(dims)]], size))
      stop("npsigtest MPI gather returned malformed array output", call. = FALSE)
    return(lapply(seq_len(size), function(j) gathered[, j]))
  }

  if (is.list(gathered)) {
    if (!identical(length(gathered), size))
      stop("npsigtest MPI gather returned malformed list output", call. = FALSE)
    return(gathered)
  }

  chunks <- as.list(gathered)
  if (!identical(length(chunks), size))
    stop("npsigtest MPI gather returned malformed atomic output", call. = FALSE)
  chunks
}

.npRmpi_npsig_parallel_boot_values_collective <- function(boot.seeds,
                                                          worker,
                                                          comm = 1L) {
  n.boot <- length(boot.seeds)
  if (n.boot < 1L)
    return(numeric(0))

  size <- mpi.comm.size(comm)
  if (size < 2L) {
    local.idx <- seq_len(n.boot)
    return(as.numeric(worker(local.idx, boot.seeds)))
  }

  rank <- mpi.comm.rank(comm)
  local.idx <- seq.int(rank + 1L, n.boot, by = size)
  local.vals <- if (length(local.idx)) {
    as.numeric(worker(local.idx, boot.seeds))
  } else {
    numeric(0)
  }

  gathered <- mpi.gather.Robj(local.vals, root = 0L, comm = comm)
  if (rank == 0L) {
    out <- numeric(n.boot)
    gathered <- .npRmpi_npsig_gather_rank_chunks(gathered = gathered, size = size)
    for (r in seq_len(size)) {
      idx.r <- seq.int(r, n.boot, by = size)
      vals.r <- as.numeric(gathered[[r]])
      if (length(idx.r) != length(vals.r))
        stop("npsigtest MPI gather returned mismatched bootstrap chunk lengths", call. = FALSE)
      if (length(idx.r))
        out[idx.r] <- vals.r
    }
    mpi.bcast.Robj(out, rank = 0L, comm = comm)
    out
  } else {
    mpi.bcast.Robj(rank = 0L, comm = comm)
  }
}

.npRmpi_npsig_bootstrap_tasks <- function(n.boot, chunk.size) {
  starts <- seq.int(1L, n.boot, by = chunk.size)
  lapply(starts, function(start) {
    list(
      start = as.integer(start),
      bsz = as.integer(min(chunk.size, n.boot - start + 1L))
    )
  })
}

.npRmpi_npsig_parallel_boot_values <- function(boot.seeds,
                                               worker,
                                               required.bindings = NULL,
                                               what = "npsigtest",
                                               profile.where = NA_character_,
                                               comm = 1L) {
  n.boot <- length(boot.seeds)
  if (n.boot < 1L)
    return(numeric(0))

  if (.npRmpi_npsig_collective_context()) {
    return(.npRmpi_npsig_parallel_boot_values_collective(
      boot.seeds = boot.seeds,
      worker = worker,
      comm = comm
    ))
  }

  if (!isTRUE(.npRmpi_has_active_slave_pool(comm = comm))) {
    return(as.numeric(worker(seq_len(n.boot), boot.seeds)))
  }

  workers <- max(1L, .npRmpi_bootstrap_worker_count(comm = comm))
  chunk.size <- max(1L, as.integer(floor(n.boot / workers)))
  chunk.size <- .npRmpi_bootstrap_tune_chunk_size(
    B = n.boot,
    chunk.size = chunk.size,
    comm = comm,
    include.master = TRUE
  )
  tasks <- .npRmpi_npsig_bootstrap_tasks(n.boot = n.boot, chunk.size = chunk.size)

  bindings <- c(list(.npsig_worker = worker), required.bindings)

  worker.chunk <- function(task, boot.seeds) {
    idx <- seq.int(task$start, length.out = task$bsz)
    vals <- as.numeric(.npsig_worker(idx, boot.seeds))
    matrix(vals, nrow = task$bsz, ncol = 1L)
  }

  out <- .npRmpi_bootstrap_run_fanout(
    tasks = tasks,
    worker = worker.chunk,
    ncol.out = 1L,
    what = what,
    profile.where = profile.where,
    comm = comm,
    master_local_chunk = TRUE,
    required.bindings = bindings,
    boot.seeds = boot.seeds
  )

  as.numeric(out[, 1L])
}

.npRmpi_npsig_bootstrap_bw_reselect <- function(xdat,
                                                ydat,
                                                bws.seed,
                                                extra.args = list(),
                                                bootstrap.iter,
                                                bw.fun = npregbw,
                                                localize = TRUE) {
  bw.args <- if (length(extra.args)) extra.args else list()
  bw.args[c("xdat", "ydat", "bws")] <- NULL

  user.nmulti <- !is.null(names(bw.args)) &&
    "nmulti" %in% names(bw.args) &&
    !is.null(bw.args$nmulti)

  if (!user.nmulti && bootstrap.iter > 1L)
    bw.args$nmulti <- 1L

  call.args <- c(list(xdat = xdat, ydat = ydat, bws = bws.seed), bw.args)

  result <- .np_progress_with_legacy_suppressed(
    if (localize) {
      .npRmpi_with_local_regression(do.call(bw.fun, call.args))
    } else {
      do.call(bw.fun, call.args)
    }
  )

  .npRmpi_autodispatch_untag(result)
}

npsigtest <-
  function(bws, ...){
    args <- list(...)

    if (!missing(bws)){
      if (is.recursive(bws)){
        if (!is.null(bws$formula) && is.null(args$xdat))
          UseMethod("npsigtest",bws$formula)
#        else if (!is.null(bws$call) && is.null(args$xdat) && (class(bws) != "npregression"))
        else if (!is.null(bws$call) && is.null(args$xdat) && (!isa(bws,"npregression")))
          UseMethod("npsigtest",bws$call)
        else if (!is.call(bws))
          UseMethod("npsigtest",bws)
        else
          UseMethod("npsigtest",NULL)
      } else {
        UseMethod("npsigtest", NULL)
      }
    } else {
      UseMethod("npsigtest", NULL)
    }
  }

npsigtest.formula <-
  function(bws, data = NULL, ...){

    tt <- terms(bws)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(bws$call), nomatch = 0)
    tmf <- bws$call[c(1,m)]
    tmf[[1]] <- as.name("model.frame")
    tmf[["formula"]] <- tt
    mf.args <- as.list(tmf)[-1L]
    umf <- tmf <- do.call(stats::model.frame, mf.args, envir = environment(tt))

    ydat <- model.response(tmf)
    xdat <- tmf[, attr(attr(tmf, "terms"),"term.labels"), drop = FALSE]

    ev <- npsigtest(xdat = xdat, ydat = ydat, bws = bws, ...)
    ev$call <- match.call(expand.dots = FALSE)
    environment(ev$call) <- parent.frame()
    ev$rows.omit <- as.vector(attr(umf,"na.action"))
    ev$nobs.omit <- length(ev$rows.omit)
    return(ev)
  }

npsigtest.call <-
  function(bws, ...) {
    ev <- npsigtest(xdat = .np_eval_bws_call_arg(bws, "xdat"),
                    ydat = .np_eval_bws_call_arg(bws, "ydat"),
                    bws = bws, ...)
    ev$call <- match.call(expand.dots = FALSE)
    environment(ev$call) <- parent.frame()
    return(ev)
  }

npsigtest.npregression <-
  function(bws, ...){
    ev <- npsigtest(bws$bws, ...)
    ev$call <- match.call(expand.dots = FALSE)
    environment(ev$call) <- parent.frame()
    return(ev)
  }

npsigtest.rbandwidth <- function(bws,
                                 xdat = stop("data xdat missing"),
                                 ydat = stop("data ydat missing"),
                                 boot.num = 399,
                                 boot.method = c("iid","wild","wild-rademacher","pairwise"),
                                 boot.type = c("I","II"),
                                 pivot = TRUE,
                                 joint = FALSE,
                                 index = seq_len(ncol(xdat)),
                                 random.seed = 42,
                                 ...) {
  .npRmpi_require_active_slave_pool(where = "npsigtest()")

  xdat <- toFrame(xdat)

  if(boot.num < 9) stop("number of bootstrap replications must be >= 9")

  ## catch and destroy NA's
  goodrows <- seq_len(nrow(xdat))
  rows.omit <- attr(na.omit(data.frame(xdat,ydat)), "na.action")
  goodrows[rows.omit] <- 0

  if (all(goodrows==0))
    stop("Data has no rows without NAs")

  xdat <- xdat[goodrows,,drop = FALSE]
  ydat <- ydat[goodrows]

  ## Fast-fail contract for invalid index values must run before any
  ## expensive/distributed bandwidth/regression work.
  .npRmpi_npsig_validate_index(index = index, xdat = xdat)

  bws <- local({
    old.disable <- getOption("npRmpi.autodispatch.disable", FALSE)
    options(npRmpi.autodispatch.disable = TRUE)
    on.exit(options(npRmpi.autodispatch.disable = old.disable), add = TRUE)
    npregbw(xdat = xdat, ydat = ydat, bws = bws, bandwidth.compute = FALSE)
  })
  bws <- .npRmpi_autodispatch_untag(bws)

  if (is.factor(ydat))
    stop("dependent variable must be continuous.")

  ## Save seed prior to setting

  seed.state <- .np_seed_enter(random.seed)


  boot.type <- match.arg(boot.type)
  boot.method <- match.arg(boot.method)
  collective.mode <- .npRmpi_npsig_collective_context()
  if (boot.type == "II")
    bws.original <- bws

  num.obs <- nrow(xdat)
  extra.args <- list(...)
  npreg.eval.fun <- if (boot.type == "II") .npRmpi_npsig_npreg_leaf else .npRmpi_npsig_do_local

  if(!joint) {

    In <- numeric(length(index))
    P <- numeric(length(index))

  }

  ## Some constants for the wild bootstrap

  a <- -0.6180339887499  # (1-sqrt(5))/2
  b <- 1.6180339887499   # (1+sqrt(5))/2
  P.a <-0.72360679774998 # (1+sqrt(5))/(2*sqrt(5))

  draw.wild.mult <- function(n.obs, a, b, p.a) {
    u <- stats::runif(n.obs)
    mult <- rep.int(b, n.obs)
    mult[u <= p.a] <- a
    mult
  }

  ## A vector for storing the resampled statistics

  In.vec <- numeric(boot.num)

  if(joint==TRUE) {

    ## Joint test

    In.mat = matrix(data = 0, ncol = 1, nrow = boot.num)

    if (boot.type == "II")
      bws <- bws.original

    ## Note - xdat must be a data frame

    ## Construct In, the average value of the squared derivatives of
    ## the jth element, discrete or continuous

    npreg.out <- npreg.eval.fun(extra.args,
                                txdat = xdat,
                                tydat = ydat,
                                bws = bws,
                                gradients = TRUE)

    In <- if(!pivot) {
      mean(npreg.out$grad[,index]^2)
    } else {
      ## Temporarily trap NaN XXX
      npreg.out$gerr[is.nan(npreg.out$gerr)] <- .Machine$double.xmax
      mean((npreg.out$grad[,index]/NZD(npreg.out$gerr[,index]))^2)
    }

    if(boot.method != "pairwise") {

      ## Compute scale and mean of unrestricted residuals

      npreg.unres <- npreg.eval.fun(extra.args,
                                    txdat = xdat,
                                    tydat = ydat,
                                    bws = bws,
                                    residuals = TRUE)
      ei.unres <- scale(npreg.unres$resid)
      ei.unres.scale <- attr(ei.unres,"scaled:scale")
      ei.unres.center <- attr(ei.unres,"scaled:center")      

      ## We now construct mhat.xi holding constant the variable whose
      ## significance is being tested at its median. First, make a copy
      ## of the data frame xdat
      
      xdat.eval <- xdat
      
      ## Impose the null by evaluating the conditional mean holding
      ## xdat[,i] constant at its median (numeric) or mode
      ## (factor/ordered) using uocquantile()

      for(i in index) {
        xq <- uocquantile(xdat[,i], 0.5)
        if (is.factor(xdat[,i]) || is.ordered(xdat[,i])) {
          xdat.eval[,i] <- cast(xq, xdat[,i], same.levels = TRUE)
        } else {
          xdat.eval[,i] <- xq
        }
      }
      
      mhat.xi <-  npreg.eval.fun(extra.args,
                                 txdat = xdat,
                                 tydat = ydat,
                                 exdat = xdat.eval,
                                 bws = bws)$mean

      ## Rescale and recenter the residuals under the null to those
      ## under the alternative
      
      ei <- as.numeric(scale(ydat-mhat.xi)*ei.unres.scale+ei.unres.center)
      
      ## Recenter the residuals to have mean zero

      ei <- ei - mean(ei)
      
    }

    .np_progress_note("Testing joint significance")
    progress <- .np_progress_begin("Bootstrap replications", total = boot.num, surface = "bootstrap")
    
    if (boot.type == "II") {
      bws.boot.prev <- bws.original

      for (i.star in seq_len(boot.num)) {
        if (boot.method == "iid") {
          ydat.star <- mhat.xi + ei[sample.int(num.obs, replace = TRUE)]
        } else if (boot.method == "wild") {
          ydat.star <- mhat.xi + ei * draw.wild.mult(num.obs, a, b, P.a)
        } else if (boot.method == "wild-rademacher") {
          ydat.star <- mhat.xi + ei * draw.wild.mult(num.obs, -1, 1, P.a)
        } else {
          boot.index <- sample.int(num.obs, replace = TRUE)
          ydat.star <- ydat[boot.index]
          xdat.star <- xdat[boot.index,]
          for (jj in index)
            xdat.star[, jj] <- xdat[, jj]
        }

        if (boot.method == "pairwise") {
          bws.boot <- .npRmpi_npsig_bootstrap_bw_reselect(
            xdat = xdat.star,
            ydat = ydat.star,
            bws.seed = bws.boot.prev,
            extra.args = extra.args,
            bootstrap.iter = i.star,
            localize = FALSE
          )
        } else {
          bws.boot <- .npRmpi_npsig_bootstrap_bw_reselect(
            xdat = xdat,
            ydat = ydat.star,
            bws.seed = bws.boot.prev,
            extra.args = extra.args,
            bootstrap.iter = i.star,
            localize = FALSE
          )
        }

        bws.boot.prev <- bws.boot
        bws <- bws.original
        bws$bw[index] <- bws.boot$bw[index]

        if (boot.method == "pairwise") {
          npreg.boot <- .npRmpi_npsig_npreg_leaf(extra.args,
                                                 txdat = xdat.star,
                                                 tydat = ydat.star,
                                                 bws = bws,
                                                 gradients = TRUE)
        } else {
          npreg.boot <- .npRmpi_npsig_npreg_leaf(extra.args,
                                                 txdat = xdat,
                                                 tydat = ydat.star,
                                                 bws = bws,
                                                 gradients = TRUE)
        }

        In.vec[i.star] <- if (!pivot) {
          mean(npreg.boot$grad[, index]^2)
        } else {
          npreg.boot$gerr[is.nan(npreg.boot$gerr)] <- .Machine$double.xmax
          mean((npreg.boot$grad[, index] / NZD(npreg.boot$gerr[, index]))^2)
        }
        progress <- .np_progress_step(progress, done = i.star)
      }
    } else {
      boot.seeds <- .npRmpi_npsig_bootstrap_seed_plan(
        num.obs = num.obs,
        boot.num = boot.num,
        boot.method = boot.method,
        draw.wild.mult = draw.wild.mult,
        a = a,
        b = b,
        p.a = P.a
      )

      joint.eval <- function(task.idx, seed.plan) {
        out <- numeric(length(task.idx))
        for (kk in seq_along(task.idx)) {
          assign(".Random.seed", seed.plan[[task.idx[kk]]], envir = .GlobalEnv)
          if (boot.method == "iid") {
            ydat.star <- mhat.xi + ei[sample.int(num.obs, replace = TRUE)]
            npreg.boot <- .npRmpi_npsig_do_local(extra.args,
                                                 txdat = xdat,
                                                 tydat = ydat.star,
                                                 bws = bws,
                                                 gradients = TRUE)
          } else if (boot.method == "wild") {
            ydat.star <- mhat.xi + ei * draw.wild.mult(num.obs, a, b, P.a)
            npreg.boot <- .npRmpi_npsig_do_local(extra.args,
                                                 txdat = xdat,
                                                 tydat = ydat.star,
                                                 bws = bws,
                                                 gradients = TRUE)
          } else if (boot.method == "wild-rademacher") {
            ydat.star <- mhat.xi + ei * draw.wild.mult(num.obs, -1, 1, P.a)
            npreg.boot <- .npRmpi_npsig_do_local(extra.args,
                                                 txdat = xdat,
                                                 tydat = ydat.star,
                                                 bws = bws,
                                                 gradients = TRUE)
          } else {
            boot.index <- sample.int(num.obs, replace = TRUE)
            ydat.star <- ydat[boot.index]
            xdat.star <- xdat[boot.index,]
            for (jj in index)
              xdat.star[, jj] <- xdat[, jj]
            npreg.boot <- .npRmpi_npsig_do_local(extra.args,
                                                 txdat = xdat.star,
                                                 tydat = ydat.star,
                                                 bws = bws,
                                                 gradients = TRUE)
          }

          out[kk] <- if (!pivot) {
            mean(npreg.boot$grad[, index]^2)
          } else {
            npreg.boot$gerr[is.nan(npreg.boot$gerr)] <- .Machine$double.xmax
            mean((npreg.boot$grad[, index] / NZD(npreg.boot$gerr[, index]))^2)
          }
        }
        out
      }

      In.vec <- .npRmpi_npsig_parallel_boot_values(
        boot.seeds = boot.seeds,
        worker = joint.eval,
        required.bindings = list(
          boot.method = boot.method,
          mhat.xi = mhat.xi,
          ei = ei,
          xdat = xdat,
          ydat = ydat,
          bws = bws,
          index = index,
          pivot = pivot,
          num.obs = num.obs,
          draw.wild.mult = draw.wild.mult,
          a = a,
          b = b,
          P.a = P.a,
          extra.args = extra.args
        ),
        what = "npsigtest",
        profile.where = "npsigtest:joint"
      )
    }

    progress <- .np_progress_end(progress)

    ## Compute the P-value

    P <- mean(In.vec > In)

    In.mat[,1] = In.vec

  } else {

    ## Individual test

    ## ii is the counter for successive elements of In and P...

    In.mat = matrix(data = 0, ncol = length(index), nrow = boot.num)

    ii <- 0

    for(i in index) {
      
      ## Increment counter...
      
      ii <- ii + 1

      if (boot.type == "II")
        bws <- bws.original
      
      ## Note - xdat must be a data frame
      
      ## Construct In, the average value of the squared derivatives of
      ## the jth element, discrete or continuous
      
      npreg.out <- npreg.eval.fun(extra.args,
                                  txdat = xdat,
                                  tydat = ydat,
                                  bws = bws,
                                  gradients = TRUE)
      
      In[ii] <- if(!pivot) {
        mean(npreg.out$grad[,i]^2)
      } else {
        ## Temporarily trap NaN XXX
        npreg.out$gerr[is.nan(npreg.out$gerr)] <- .Machine$double.xmax
        mean((npreg.out$grad[,i]/NZD(npreg.out$gerr[,i]))^2)
      }
      
      if(boot.method != "pairwise") {

        ## Compute scale and mean of unrestricted residuals

        npreg.unres <- npreg.eval.fun(extra.args,
                                      txdat = xdat,
                                      tydat = ydat,
                                      bws = bws,
                                      residuals = TRUE)
        ei.unres <- scale(npreg.unres$resid)
        ei.unres.scale <- attr(ei.unres,"scaled:scale")
        ei.unres.center <- attr(ei.unres,"scaled:center")      

        ## We now construct mhat.xi holding constant the variable whose
        ## significance is being tested at its median. First, make a copy
        ## of the data frame xdat
        
        xdat.eval <- xdat
        
        ## Impose the null by evaluating the conditional mean holding
        ## xdat[,i] constant at its median (numeric) or mode
        ## (factor/ordered) using uocquantile()
        
        xq <- uocquantile(xdat[,i], 0.5)
        if (is.factor(xdat[,i]) || is.ordered(xdat[,i])) {
          xdat.eval[,i] <- cast(xq, xdat[,i], same.levels = TRUE)
        } else {
          xdat.eval[,i] <- xq
        }
        
        mhat.xi <-  npreg.eval.fun(extra.args,
                                   txdat = xdat,
                                   tydat = ydat,
                                   exdat = xdat.eval,
                                   bws = bws)$mean
        
        ## Rescale and recenter the residuals under the null to those
        ## under the alternative
        
        ei <- as.numeric(scale(ydat-mhat.xi)*ei.unres.scale+ei.unres.center)
        
        ## Recenter the residuals to have mean zero
        
        ei <- ei - mean(ei)
        
      }

      .np_progress_note(sprintf("Testing variable %s of (%s)", i, paste(index, collapse = ",")))
      progress <- .np_progress_begin("Bootstrap replications", total = boot.num, surface = "bootstrap")
      
      if (boot.type == "II") {
        bws.boot.prev <- bws.original

        for (i.star in seq_len(boot.num)) {
          if (boot.method == "iid") {
            ydat.star <- mhat.xi + ei[sample.int(num.obs, replace = TRUE)]
          } else if (boot.method == "wild") {
            ydat.star <- mhat.xi + ei * draw.wild.mult(num.obs, a, b, P.a)
          } else if (boot.method == "wild-rademacher") {
            ydat.star <- mhat.xi + ei * draw.wild.mult(num.obs, -1, 1, P.a)
          } else {
            boot.index <- sample.int(num.obs, replace = TRUE)
            ydat.star <- ydat[boot.index]
            xdat.star <- xdat
            xdat.star[, -i] <- xdat[boot.index, -i]
          }

          if (boot.method == "pairwise") {
            bws.boot <- .npRmpi_npsig_bootstrap_bw_reselect(
              xdat = xdat.star,
              ydat = ydat.star,
              bws.seed = bws.boot.prev,
              extra.args = extra.args,
              bootstrap.iter = i.star,
              localize = FALSE
            )
          } else {
            bws.boot <- .npRmpi_npsig_bootstrap_bw_reselect(
              xdat = xdat,
              ydat = ydat.star,
              bws.seed = bws.boot.prev,
              extra.args = extra.args,
              bootstrap.iter = i.star,
              localize = FALSE
            )
          }

          bws.boot.prev <- bws.boot
          bws <- bws.original
          bws$bw[i] <- bws.boot$bw[i]

          if (boot.method == "pairwise") {
            npreg.boot <- .npRmpi_npsig_npreg_leaf(extra.args,
                                                   txdat = xdat.star,
                                                   tydat = ydat.star,
                                                   bws = bws,
                                                   gradients = TRUE)
          } else {
            npreg.boot <- .npRmpi_npsig_npreg_leaf(extra.args,
                                                   txdat = xdat,
                                                   tydat = ydat.star,
                                                   bws = bws,
                                                   gradients = TRUE)
          }

          In.vec[i.star] <- if (!pivot) {
            mean(npreg.boot$grad[, i]^2)
          } else {
            npreg.boot$gerr[is.nan(npreg.boot$gerr)] <- .Machine$double.xmax
            mean((npreg.boot$grad[, i] / NZD(npreg.boot$gerr[, i]))^2)
          }
          progress <- .np_progress_step(progress, done = i.star)
        }
      } else {
        boot.seeds <- .npRmpi_npsig_bootstrap_seed_plan(
          num.obs = num.obs,
          boot.num = boot.num,
          boot.method = boot.method,
          draw.wild.mult = draw.wild.mult,
          a = a,
          b = b,
          p.a = P.a
        )

        indiv.eval <- function(task.idx, seed.plan) {
          out <- numeric(length(task.idx))
          for (kk in seq_along(task.idx)) {
            assign(".Random.seed", seed.plan[[task.idx[kk]]], envir = .GlobalEnv)
            if (boot.method == "iid") {
              ydat.star <- mhat.xi + ei[sample.int(num.obs, replace = TRUE)]
              npreg.boot <- .npRmpi_npsig_do_local(extra.args,
                                                   txdat = xdat,
                                                   tydat = ydat.star,
                                                   bws = bws,
                                                   gradients = TRUE)
            } else if (boot.method == "wild") {
              ydat.star <- mhat.xi + ei * draw.wild.mult(num.obs, a, b, P.a)
              npreg.boot <- .npRmpi_npsig_do_local(extra.args,
                                                   txdat = xdat,
                                                   tydat = ydat.star,
                                                   bws = bws,
                                                   gradients = TRUE)
            } else if (boot.method == "wild-rademacher") {
              ydat.star <- mhat.xi + ei * draw.wild.mult(num.obs, -1, 1, P.a)
              npreg.boot <- .npRmpi_npsig_do_local(extra.args,
                                                   txdat = xdat,
                                                   tydat = ydat.star,
                                                   bws = bws,
                                                   gradients = TRUE)
            } else {
              boot.index <- sample.int(num.obs, replace = TRUE)
              ydat.star <- ydat[boot.index]
              xdat.star <- xdat
              xdat.star[, -i] <- xdat[boot.index, -i]
              npreg.boot <- .npRmpi_npsig_do_local(extra.args,
                                                   txdat = xdat.star,
                                                   tydat = ydat.star,
                                                   bws = bws,
                                                   gradients = TRUE)
            }

            out[kk] <- if (!pivot) {
              mean(npreg.boot$grad[, i]^2)
            } else {
              npreg.boot$gerr[is.nan(npreg.boot$gerr)] <- .Machine$double.xmax
              mean((npreg.boot$grad[, i] / NZD(npreg.boot$gerr[, i]))^2)
            }
          }
          out
        }

        In.vec <- .npRmpi_npsig_parallel_boot_values(
          boot.seeds = boot.seeds,
          worker = indiv.eval,
          required.bindings = list(
            boot.method = boot.method,
            mhat.xi = mhat.xi,
            ei = ei,
            xdat = xdat,
            ydat = ydat,
            bws = bws,
            i = i,
            pivot = pivot,
            num.obs = num.obs,
            draw.wild.mult = draw.wild.mult,
            a = a,
            b = b,
            P.a = P.a,
            extra.args = extra.args
          ),
          what = "npsigtest",
          profile.where = "npsigtest:indiv"
        )
      }

      progress <- .np_progress_end(progress)
      
      ## Compute the P-value
      
      P[ii] <- mean(In.vec > In[ii])
      
      In.mat[,ii] = In.vec
      
    }
    
  } ## End invididual test

  ## Return a list containing the statistic and its P-value
  ## bootstrapped In.vec for each variable...

  ## Restore seed

  .np_seed_exit(seed.state)

  sigtest(In=In,
          In.mat,
          P=P,
          bws = bws,
          ixvar = index,
          boot.method,
          pivot,
          joint,
          boot.type,
          boot.num)

}

npsigtest.default <- function(bws, xdat, ydat, ...){
  .npRmpi_require_active_slave_pool(where = "npsigtest()")

  sc <- sys.call()
  sc.names <- names(sc)

  ## here we check to see if the function was called with tdat = if it
  ## was, we need to catch that and map it to dat = otherwise the call
  ## is passed unadulterated to npudensbw

  bws.named <- any(sc.names == "bws")
  xdat.named <- any(sc.names == "xdat")
  ydat.named <- any(sc.names == "ydat")

  no.bws <- missing(bws)
  no.xdat <- missing(xdat)
  no.ydat <- missing(ydat)

  ## if bws was passed in explicitly, do not compute bandwidths

  ## autodispatch normalizes calls via match.call(), which can turn an
  ## originally unnamed formula first argument into named bws=... .
  ## Preserve legacy formula behavior by rewriting npregbw() call shape.
  if (bws.named && no.xdat && no.ydat && inherits(bws, "formula")) {
    sc$`bws` <- NULL
    sc$formula <- bws
    sc.bw <- sc
    sc.bw[[1]] <- quote(npregbw)
    bws.named <- FALSE
  } else {
    sc.bw <- sc
    sc.bw[[1]] <- quote(npregbw)
  }

  if(xdat.named)
    xdat <- toFrame(xdat)

  if(bws.named){
    sc.bw$bandwidth.compute <- FALSE
  }

  tbw <- .np_eval_bw_call(sc.bw, caller_env = parent.frame())
  
  call.args <- list(bws = tbw)
  if(!no.xdat)
    call.args$xdat <- xdat
  if(!no.ydat)
    call.args$ydat <- ydat

  dots <- list(...)
  dots[c("bws", "bandwidth.compute", "formula", "data", "xdat", "ydat")] <- NULL

  ev <- do.call(npsigtest, c(call.args, dots))

  ev$call <- match.call(expand.dots = FALSE)
  environment(ev$call) <- parent.frame()
  return(ev)
}
