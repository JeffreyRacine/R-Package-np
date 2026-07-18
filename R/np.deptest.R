## Function that implements the bivariate dependence metric described
## in Maasoumi, E.~and J.S.~Racine (2002), "Entropy and Predictability
## of Stock Market Returns," Journal of Econometrics, March, Volume
## 107, Issue 2, pp 291-312.

.npRmpi_dept_collective_context <- function() {
  isTRUE(.npRmpi_autodispatch_called_from_bcast())
}

.npRmpi_dept_bootstrap_index_plan <- function(n, boot.num) {
  index <- matrix(NA_integer_, nrow = boot.num, ncol = n)

  for (b in seq_len(boot.num))
    index[b, ] <- sample.int(n, replace = TRUE)

  index
}

.npRmpi_dept_numeric_chunks <- function(gathered, size) {
  size <- as.integer(size)[1L]
  if (is.na(size) || size < 1L)
    stop("invalid MPI gather size")

  if (is.matrix(gathered)) {
    if (!identical(ncol(gathered), size))
      stop("npdeptest MPI gather returned malformed matrix output", call. = FALSE)
    return(lapply(seq_len(size), function(j) gathered[, j]))
  }

  if (is.array(gathered)) {
    dims <- dim(gathered)
    if (length(dims) < 2L || !identical(dims[[length(dims)]], size))
      stop("npdeptest MPI gather returned malformed array output", call. = FALSE)
    return(lapply(seq_len(size), function(j) gathered[, j]))
  }

  if (is.list(gathered)) {
    if (!identical(length(gathered), size))
      stop("npdeptest MPI gather returned malformed list output", call. = FALSE)
    return(gathered)
  }

  chunks <- as.list(gathered)
  if (!identical(length(chunks), size))
    stop("npdeptest MPI gather returned malformed atomic output", call. = FALSE)
  chunks
}

.npRmpi_dept_collective_bootstrap <- function(plan,
                                              data.x,
                                              data.y,
                                              bw.data.x,
                                              bw.data.y,
                                              bw.joint,
                                              Srho.bivar,
                                              method,
                                              progress = NULL,
                                              comm = 1L) {
  boot.num <- nrow(plan)
  size <- mpi.comm.size(comm)
  rank <- mpi.comm.rank(comm)

  local.idx <- seq.int(rank + 1L, boot.num, by = size)
  .npRmpi_bootstrap_transport_trace(
    what = "npdeptest",
    event = "fanout.collective.start",
    fields = list(rank = rank, size = size, B = boot.num, local = length(local.idx))
  )

  local.Srho <- numeric(length(local.idx))

  for (jj in seq_along(local.idx)) {
    b <- local.idx[[jj]]
    data.x.boot <- data.x[plan[b, ]]
    local.Srho[[jj]] <- .npRmpi_with_local_regression(
      Srho.bivar(data.x.boot, data.y, bw.data.x, bw.data.y, bw.joint, method = method)
    )
  }

  invisible(gc(FALSE))

  payload <- c(as.numeric(local.idx), local.Srho)
  gathered <- mpi.gather.Robj(payload, root = 0L, comm = comm)

  if (rank == 0L) {
    chunks <- .npRmpi_dept_numeric_chunks(gathered, size = size)
    Srho.vec.boot <- numeric(boot.num)

    for (rr in seq_len(size)) {
      vals <- as.numeric(chunks[[rr]])
      if (!length(vals))
        next
      if ((length(vals) %% 2L) != 0L)
        stop("npdeptest MPI gather returned malformed bootstrap chunk", call. = FALSE)
      n.local <- length(vals) / 2L
      idx <- as.integer(vals[seq_len(n.local)])
      Srho.vec.boot[idx] <- vals[n.local + seq_len(n.local)]
    }

    if (!is.null(progress))
      progress <- .np_progress_step(progress, done = boot.num)
    .npRmpi_bootstrap_transport_trace(
      what = "npdeptest",
      event = "fanout.collective.done",
      fields = list(rank = rank, size = size, B = boot.num)
    )
    mpi.bcast.Robj(Srho.vec.boot, rank = 0L, comm = comm)
    Srho.vec.boot
  } else {
    .npRmpi_bootstrap_transport_trace(
      what = "npdeptest",
      event = "fanout.collective.done",
      fields = list(rank = rank, size = size, B = boot.num)
    )
    mpi.bcast.Robj(rank = 0L, comm = comm)
  }
}

npdeptest <- function(data.x = NULL,
                      data.y = NULL,
                      method=c("integration","summation"),
                      bootstrap = TRUE,
                      boot.num = 399,
                      random.seed = 42) {
  .npRmpi_require_active_slave_pool(where = "npdeptest()")
  if (.npRmpi_autodispatch_active())
    return(.npRmpi_autodispatch_call(match.call(), parent.frame()))
  
  ## Trap fatal errors

  if(is.data.frame(data.x)||is.data.frame(data.y)) stop(" you must enter two data vectors (and not data frames)")
  if(is.factor(data.x)||is.factor(data.y)) stop(" does not support factors")
  if(is.null(data.x)||is.null(data.y)) stop(" you must enter x and y data vectors")
  if(ncol(data.frame(data.x)) != 1) stop(" data must have one dimension only")
  if(length(data.x)!=length(data.y)) stop(" data vectors must be of equal length")
  if(boot.num < 9) stop(" number of bootstrap replications must be >= 9")

  method <- match.arg(method)

  ## Save seed prior to setting

  seed.state <- .np_seed_enter(random.seed)
  on.exit(.np_seed_exit(seed.state, remove_if_absent = TRUE), add = TRUE)

  .np_progress_note("Computing bandwidths")

  ## If the variable is a time series convert to type numeric

  if(is.ts(data.x)) data.x <- as.numeric(data.x)
  if(is.ts(data.y)) data.y <- as.numeric(data.y)

  ## Remove any NAs from paired data

  tmp <- na.omit(data.frame(data.x,data.y))
  data.x <- tmp$data.x
  data.y <- tmp$data.y
  rm(tmp)

  ## Define the metric entropy function

  Srho.bivar <- function(x.dat = NULL,
                         y.dat = NULL,
                         bw.x = NULL,
                         bw.y = NULL,
                         bw.joint = NULL,
                         method=c("integration","summation")) {
    
    if(is.null(bw.x)||is.null(bw.y)||is.null(bw.joint)) stop(" you must provide numeric bandwidths for f(x), f(y) and f(x,y)")
    
    method <- match.arg(method)
    
    if(method=="summation") {
      
      ## Summation version: \sum_i(1-sqrt(f(x_i)f(y_i)/f(x_i,y_i)))^2
      ## We could code this by hand per below, however, there is no
      ## performance penalty imposed by sum() so no need.
      
      f.x <- fitted(npudens(tdat=x.dat,bws=bw.x))
      f.y <- fitted(npudens(tdat=y.dat,bws=bw.y))
      f.xy <- fitted(npudens(tdat=cbind(x.dat,y.dat),bws=bw.joint))
      summand <- f.x*f.y/f.xy
      
      ## In summation version we divide densities which can lead to
      ## numerical instability issues not present in the integrand
      ## version (which uses differences instead). We check for this
      ## case, remove offending points, and warn. This traps -Inf,
      ## Inf, and NaN.  

      if(!all(is.finite(summand))) {
        .np_warning(" non-finite value in summation-based statistic: integration recommended")
        summand <- summand[is.finite(summand)]
      }

      return(0.5*mean((1-sqrt(summand))**2))
      
    } else {
      
      ## Integration version:
      ## \int\int (sqrt(f(x,y))-sqrt(f(x))sqrt(f(y)))^2 dx dy
      
      lo.default <- c(min(x.dat)-10.0*IQR(x.dat),min(y.dat)-10.0*IQR(y.dat))
      up.default <- c(max(x.dat)+10.0*IQR(x.dat),max(y.dat)+10.0*IQR(y.dat))
      
      return(.np_entropy_bivariate_integral(
        x.dat = x.dat,
        y.dat = y.dat,
        bw.x = bw.x,
        bw.y = bw.y,
        bw.joint = bw.joint,
        lower = lo.default,
        upper = up.default
      ))
    }

  } ## end of Srho.bivar function
  
  ## Compute and save bandwidths (save for bootstrapping if requested)

  bw.data.x <- .np_progress_with_legacy_suppressed(npudensbw(~data.x))$bw
  bw.data.y <- .np_progress_with_legacy_suppressed(npudensbw(~data.y))$bw
  bw.joint <- .np_progress_with_legacy_suppressed(npudensbw(~data.x+data.y))$bw

  .np_progress_note("Constructing metric entropy")
  
  Srho.vec <- Srho.bivar(data.x,data.y,bw.data.x,bw.data.y,bw.joint,method=method)

  ## Bootstrap if requested - null is independence so simple iid
  ## index resampling under replacement is sufficient

  if(bootstrap) {

    progress <- .np_progress_begin("Bootstrap replications", total = boot.num, surface = "bootstrap")

    if (.npRmpi_dept_collective_context()) {
      plan <- .npRmpi_dept_bootstrap_index_plan(length(data.x), boot.num)
      post.boot.seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      Srho.vec.boot <- .npRmpi_dept_collective_bootstrap(
        plan = plan,
        data.x = data.x,
        data.y = data.y,
        bw.data.x = bw.data.x,
        bw.data.y = bw.data.y,
        bw.joint = bw.joint,
        Srho.bivar = Srho.bivar,
        method = method,
        progress = progress
      )
      assign(".Random.seed", post.boot.seed, envir = .GlobalEnv)
    } else {
      Srho.vec.boot <- numeric()
      for (b in seq_len(boot.num)) {
        ## Break systematic relationship between x and y (null)

        data.x.boot <- data.x[sample.int(length(data.x), replace = TRUE)]

        Srho.vec.boot[b] <- Srho.bivar(data.x.boot,data.y,bw.data.x,bw.data.y,bw.joint,method=method)
        progress <- .np_progress_step(progress, done = b)
      }
    }

    progress <- .np_progress_end(progress)

    ## Compute P-values

    P <- mean(Srho.vec.boot > Srho.vec)

  }

  ## Restore seed

  .np_seed_exit(seed.state, remove_if_absent = TRUE)

  if(bootstrap) {

    deptest(Srho = Srho.vec,
            Srho.bootstrap.vec = Srho.vec.boot,
            P = P,
            bootstrap = bootstrap,
            boot.num = boot.num,
            bw.data.x = bw.data.x,
            bw.data.y = bw.data.y,
            bw.joint = bw.joint)

  } else {

    deptest(Srho = Srho.vec,
            Srho.bootstrap.vec = NULL,
            P = NULL,
            bootstrap = bootstrap,
            boot.num = NULL,
            bw.data.x = bw.data.x,
            bw.data.y = bw.data.y,
            bw.joint = bw.joint)

  }

}
