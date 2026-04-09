### Copyright (C) 2002 Hao Yu 

.npRmpi_tail_slave_indices <- function(size, caller = "slave.hostinfo") {
    size <- as.integer(size)
    if (length(size) != 1L || is.na(size) || size < 2L)
        stop(sprintf("%s: 'size' must be an integer >= 2", caller))
    seq.int(max(1L, size - 2L), size - 1L)
}

mpi.hostinfo <- function(comm=1){
    if (mpi.comm.size(comm)==0){
        err <-paste("It seems no members running on comm", comm)
        stop(err)
    }
    hostname <- mpi.get.processor.name() 
    rk <- mpi.comm.rank(comm=comm)
    size <- mpi.comm.size(comm=comm)
    cat("\tHost:",hostname,"\tRank(ID):",rk, "\tof Size:", size,
        "on comm", comm, "\n")
}

slave.hostinfo <- function(comm=1, short=TRUE){
    #if (!mpi.is.master())
 	if (mpi.comm.rank(comm)!=0)
		stop("cannot run slavehostinfo on slaves")
    size <- mpi.comm.size(comm)
    if (size==0){
        err <-paste("It seems no slaves running on comm", comm)
        stop(err)
    }
    if (size == 1)
    mpi.hostinfo(comm)
    else { 
        master <-mpi.get.processor.name() 
        slavehost <- unlist(mpi.remote.exec(mpi.get.processor.name(),comm=comm))
        slavecomm <- 1 #as.integer(mpi.remote.exec(.comm,comm=comm))
        ranks <- seq_len(size - 1L)
        commm <- paste(comm, ")",sep="")
        if (size > 10){
        rank0 <- paste("master  (rank 0 , comm", commm)
            ranks <- c(
              paste(ranks[seq_len(9L)], " ", sep = ""),
              ranks[seq.int(10L, size - 1L)]
            )
        }
        else
        rank0 <- paste("master (rank 0, comm", commm)
        cat(rank0, "of size", size, "is running on:",master, "\n")
        slavename <- paste("slave", ranks,sep="")
        ranks <- paste("(rank ",ranks, ", comm ",slavecomm,")", sep="")
			if (short && size > 8){
	          for (i in seq_len(3L)) {
	            cat(slavename[i], ranks[i], "of size",size, 
	          "is running on:",slavehost[i], "\n")	
			  }
			  cat("... ... ...\n")
			  for (i in .npRmpi_tail_slave_indices(size)){
			    cat(slavename[i], ranks[i], "of size",size, 
	          "is running on:",slavehost[i], "\n")
			  }
			}
		else {
          for (i in seq_len(size - 1L)){
            cat(slavename[i], ranks[i], "of size",size, 
          "is running on:",slavehost[i], "\n")
          }
		}
    }
}

lamhosts <- function(){
    hosts <- system("lamnodes C -c -n", TRUE)
    base <-character(0)
    for (host in hosts)
        base <- c(base, unlist(strsplit(host, "\\."))[1])
    nn <- seq_len(length(hosts)) - 1L
        names(nn) <- base
    nn
}

mpi.spawn.Rslaves <- 
    function(Rscript=system.file("slavedaemon.R", package="npRmpi"),
    nslaves=mpi.universe.size(),
    root=0,
    intercomm=2,
    comm=1,
    hosts=NULL,
    needlog=FALSE,
    mapdrive=TRUE,
	quiet=FALSE,
	nonblock=TRUE,
	sleep=0.1) {
    if (!is.loaded("mpi_comm_spawn"))
        stop("You cannot use MPI_Comm_spawn API")   
    if (mpi.comm.size(comm) > 1){
         if (isTRUE(getOption("npRmpi.reuse.slaves", FALSE))) {
             existing <- mpi.comm.size(comm) - 1L
             if (!missing(nslaves) && length(nslaves) == 1 && is.finite(nslaves)) {
                 requested <- as.integer(nslaves)
                 if (!is.na(requested) && requested != existing) {
                     .np_warning(paste0(
                         "npRmpi is reusing an existing slave pool (", existing,
                         " slave(s)) but nslaves=", requested, " was requested. ",
                         "To change the number of slaves, call mpi.close.Rslaves(force=TRUE) ",
                         "and then mpi.spawn.Rslaves(nslaves=...)."
                     ), call. = FALSE)
                 }
             }
             if (!quiet) mpi.hostinfo(comm)
             return(invisible(existing))
         }
         err <-paste("It seems there are some slaves running on comm ", comm)
         stop(err)
    }
	if (.Platform$OS=="windows"){
		stop("Spawning is not implemented on Windows. Launch with mpiexec and use npRmpi.init(mode=\"attach\").")
    } else {
        tmp <- paste(Sys.getpid(), "+", comm, sep="")
        lib.path.arg <- paste(.libPaths(), collapse = .Platform$path.sep)
        if (needlog)
            arg <- c(Rscript, tmp, "needlog", R.home(), lib.path.arg)
        else
            arg <- c(Rscript, tmp , "nolog", R.home(), lib.path.arg)
        if (!is.null(hosts)){
            hosts <- as.integer(hosts)
            if (any(is.na(hosts)))
                stop("hosts argument contains non-integer object(s).")
            if (max(hosts) > mpi.universe.size() -1 ||min(hosts) < 0){
                tmp1 <- paste("hosts number should be within 0 to",
                    mpi.universe.size()-1)
                stop(tmp1)
            }
            nslaves <- length(hosts)
            tmpfile <-paste(tmp, "appschema", sep="") 
            fileobj <- file(tmpfile,"w")
            cat("c", paste(hosts, collapse=","), sep="", file=fileobj)
            cat(" ", system.file("Rslaves.sh", package="npRmpi"), file=fileobj)
            cat(" ", paste(arg, collapse=" "), file=fileobj)
            close(fileobj)
            mpi.info.create(0)
            mpi.info.set(0,"file",tmpfile)
        }
		if (length(unlist(strsplit(.Platform$pkgType,"mac"))) ==2 && .Platform$r_arch =="x86_64")
			realns<-mpi.comm.spawn(slave=system.file("MacR64slaves.sh", package="npRmpi"),
			slavearg=arg, nslaves=nslaves, info=0, root=root, intercomm=intercomm, quiet = quiet)
		else 
			realns<-mpi.comm.spawn(slave=system.file("Rslaves.sh", package="npRmpi"),
			slavearg=arg, nslaves=nslaves, info=0, root=root, intercomm=intercomm, quiet = quiet)
	}
    if (!is.null(hosts)){
        unlink(tmpfile)
        mpi.info.free(0)
    }
    if (realns==0)
        stop("It seems no single slave spawned.")
    if (mpi.intercomm.merge(intercomm,0,comm)) {
        mpi.comm.set.errhandler(comm)
        mpi.comm.disconnect(intercomm)
		mpi.bcast(nonblock,type=1, rank=0, comm=comm)
		mpi.bcast(sleep,type=2, rank=0, comm=comm)
        if (!quiet) mpi.hostinfo(comm)
    }
    else
        stop("Fail to merge the comm for master and slaves.")
}   

mpi.remote.exec <- function(cmd, ...,  simplify=TRUE, comm=1, ret=TRUE){
    if (mpi.comm.size(comm) < 2)
    stop("It seems no slaves running.")
    tag <- floor(runif(1,20000,30000))
    scmd <- substitute(cmd)
    arg <-list(...)
    .npRmpi_transport_trace(
        role = "master",
        event = "remote.exec.start",
        fields = list(tag = tag, ret = ret, simplify = simplify, comm = comm)
    )
    #if (length(arg) > 0) 
    #    deparse(arg)
    #tag.ret <- c(tag, ret, simplify)
    mpi.bcast.cmd(.mpi.worker.exec, tag=tag, ret=ret, simplify=simplify, comm = comm)
    #mpi.bcast(as.integer(tag.ret), type=1, comm=comm)
    mpi.bcast.Robj(list(scmd=scmd, arg=arg), comm=comm)

    if (ret){
        size <- mpi.comm.size(comm) 
        allcode <- mpi.allgather(integer(2), 1, integer(2*size), comm)
    .npRmpi_transport_trace(
        role = "master",
        event = "remote.exec.typecodes",
        fields = list(tag = tag, size = size)
    )
    type <- allcode[seq(3,2*size,2)]
    len <- allcode[seq(4,2*size,2)]
    eqlen <- all(len==len[1])
    if (all(type==1)){
        .npRmpi_transport_trace(
            role = "master",
            event = "remote.exec.path.int",
            fields = list(tag = tag, eqlen = eqlen, len = if (length(len)) len[1] else NA_integer_)
        )
        if (eqlen && simplify){
            out <- mpi.gather(integer(len[1]),1,integer(size*len[1]),0,comm)
            out <- out[(len[1]+1):(size*len[1])]
        dim(out) <- c(len[1], size-1)
        out <- data.frame(out)
        }
        else {
          out1<-mpi.gatherv(integer(1),1,integer(1+sum(len)),c(1,len),0,comm)
            uplen <- cumsum(len)+1
            lowlen <-c(2, uplen[-(size-1)]+1)
                out <- as.list(integer(size-1))
                names(out) <- paste("slave", seq_len(size - 1L), sep="")
                for (i in seq_len(size - 1L))
            out[[i]]<- out1[lowlen[i]:uplen[i]]
        }
    }
    else if (all(type==2)){
        .npRmpi_transport_trace(
            role = "master",
            event = "remote.exec.path.dbl",
            fields = list(tag = tag, eqlen = eqlen, len = if (length(len)) len[1] else NA_integer_)
        )
        if (eqlen && simplify){
                out <- mpi.gather(double(len[1]),2,double(size*len[1]),0,comm)
                out <- out[(len[1]+1):(size*len[1])]
                dim(out) <- c(len[1], size-1)
                out <- data.frame(out)
            }
        else {
          out1<-mpi.gatherv(double(1),2,double(1+sum(len)),c(1,len),0,comm) 
            uplen <- cumsum(len)+1
            lowlen <-c(2, uplen[-(size-1)]+1)
                out <- as.list(integer(size-1))
                names(out) <- paste("slave", seq_len(size - 1L), sep="")
                for (i in seq_len(size - 1L))
            out[[i]]<- out1[lowlen[i]:uplen[i]]
        }
    }
 else if (all(type==4)){
        .npRmpi_transport_trace(
            role = "master",
            event = "remote.exec.path.raw",
            fields = list(tag = tag, eqlen = eqlen, len = if (length(len)) len[1] else NA_integer_)
        )
        if (eqlen && simplify){
                out <- mpi.gather(raw(len[1]),4,raw(size*len[1]),0,comm)
                out <- out[(len[1]+1):(size*len[1])]
                dim(out) <- c(len[1], size-1)
                out <- data.frame(out)
            }
        else {
          out1<-mpi.gatherv(raw(1),4,raw(1+sum(len)),c(1,len),0,comm)
            uplen <- cumsum(len)+1
            lowlen <-c(2, uplen[-(size-1)]+1)
                out <- as.list(integer(size-1))
                names(out) <- paste("slave", seq_len(size - 1L), sep="")
                for (i in seq_len(size - 1L))
                    out[[i]]<- out1[lowlen[i]:uplen[i]]
        }
    }

    else {
            .npRmpi_transport_trace(
                role = "master",
                event = "remote.exec.path.robj",
                fields = list(tag = tag)
            )
            out <- as.list(integer(size-1))
            names(out) <- paste("slave", seq_len(size - 1L), sep="")
            for (i in seq_len(size - 1L)){
        tmp<- mpi.recv.Robj(mpi.any.source(),tag,comm)
        src <- mpi.get.sourcetag()[1] 
        out[[src]]<- tmp 
        }
        }
        .npRmpi_transport_trace(
            role = "master",
            event = "remote.exec.done",
            fields = list(tag = tag, ret = ret)
        )
        out
    } else {
        .npRmpi_transport_trace(
            role = "master",
            event = "remote.exec.done",
            fields = list(tag = tag, ret = ret)
        )
    }
}

.typeindex <- function (x) {
    if(inherits(x, "integer"))
            as.integer(c(1,length(x)))
    else if (inherits(x, "numeric"))
            as.integer(c(2,length(x)))
    else if (inherits(x, "raw"))
            as.integer(c(4,length(x)))

    else
            # Keep a fixed-length control packet across ranks; a scalar -1 can
            # truncate MPI_Allgather when another rank reports c(type, len).
            as.integer(c(-1L, 0L))
}

.npRmpi_transport_trace_path <- function() {
    path.opt <- getOption("npRmpi.transport.trace.file", "")
    path.env <- Sys.getenv("NP_RMPI_TRANSPORT_TRACE_FILE", unset = "")
    path <- if (nzchar(path.env)) path.env else path.opt
    path <- as.character(path)[1L]
    if (is.na(path) || !nzchar(path))
        return("")
    path
}

.npRmpi_transport_trace <- function(role, event, fields = list()) {
    path <- .npRmpi_transport_trace_path()
    if (!nzchar(path))
        return(invisible(FALSE))

    dirpath <- dirname(path)
    if (!identical(dirpath, ".") && !dir.exists(dirpath)) {
        ok <- tryCatch({
            dir.create(dirpath, recursive = TRUE, showWarnings = FALSE)
        }, error = function(e) FALSE)
        if (!isTRUE(ok) && !dir.exists(dirpath))
            return(invisible(FALSE))
    }

    if (is.null(names(fields)))
        names(fields) <- paste0("f", seq_along(fields))
    if (length(fields) > 0) {
        fields <- fields[!is.na(names(fields)) & nzchar(names(fields))]
    }
    kv <- if (length(fields) > 0) {
        paste(
            paste0(
                names(fields),
                "=",
                vapply(fields, function(v) {
                    vv <- as.character(v)[1L]
                    if (is.na(vv)) "NA" else vv
                }, character(1))
            ),
            collapse = "\t"
        )
    } else ""

    line <- paste(
        format(Sys.time(), "%Y-%m-%dT%H:%M:%OS6%z"),
        paste0("pid=", Sys.getpid()),
        paste0("role=", as.character(role)[1L]),
        paste0("event=", as.character(event)[1L]),
        kv,
        sep = "\t"
    )

    tryCatch({
        cat(paste0(line, "\n"), file = path, append = TRUE)
        TRUE
    }, error = function(e) FALSE)
}

.npRmpi_eval_scmd <- function(scmd, arg = list(), envir = parent.frame()) {
    if (length(arg) > 0)
        return(do.call(.npRmpi_bcast_cmd_funref(scmd), arg, envir = envir))
    if (is.symbol(scmd) && is.environment(envir)) {
        not_found <- new.env(parent = emptyenv())
        sym_val <- get0(as.character(scmd), envir = envir, inherits = TRUE, ifnotfound = not_found)
        if (!identical(sym_val, not_found))
            return(sym_val)
    }
    out <- .np_try_eval_in_frames(scmd, eval_env = envir, search_frames = FALSE)
    if (isTRUE(out$ok))
        return(out$value)
    if (inherits(out$error, "error"))
        stop(conditionMessage(out$error), call. = FALSE)
    stop("unable to evaluate MPI command expression", call. = FALSE)
}

.npRmpi_clear_applylb_cache <- function() {
    if (!is.null(get0(".mpi.applyLB", envir = .GlobalEnv, inherits = FALSE)))
        rm(".mpi.applyLB", envir = .GlobalEnv)
    invisible(NULL)
}

.mpi.worker.exec <- function(tag, ret, simplify){
	.comm <- 1
    .npRmpi_transport_trace(
        role = "worker",
        event = "worker.exec.start",
        fields = list(tag = tag, ret = ret, simplify = simplify, comm = .comm)
    )
    #tag.ret <- mpi.bcast(integer(3), type=1, comm=.comm)
    #tag <- tag.ret[1]
    #ret <- as.logical(tag.ret[2])
    #simplify <- as.logical(tag.ret[3])
    scmd.arg <- mpi.bcast.Robj(comm=.comm)

	    if (ret){
	    size <- mpi.comm.size(.comm)
	    out <- tryCatch(.npRmpi_eval_scmd(scmd.arg$scmd, scmd.arg$arg, envir = sys.parent()),
	                    error = function(e) e)

	    type <- .typeindex(out)
	    if (is.na(type[2]))
	        type[2] <- as.integer(0)
        .npRmpi_transport_trace(
            role = "worker",
            event = "worker.exec.type",
            fields = list(tag = tag, type = type[1], len = type[2], size = size)
        )
	    allcode <- mpi.allgather(type, 1, integer(2*size), .comm)
    type <- allcode[seq(3,2*size,2)]
        len <- allcode[seq(4,2*size,2)]
        eqlen <- all(len==len[1])
        if (all(type==1)) {
            if (eqlen && simplify)
                mpi.gather(out, 1, integer(1), 0, .comm)
        else
        mpi.gatherv(out, 1, integer(1), integer(1), 0 ,.comm)
    }
    else if (all(type==2)) {
            if (eqlen && simplify)
                mpi.gather(out, 2, double(1), 0, .comm)
        else
                mpi.gatherv(out, 2, double(1), integer(1), 0, .comm)
        }
     else if (all(type==4)) {
            if (eqlen && simplify)
                mpi.gather(out, 4, raw(1), 0, .comm)
        else
                mpi.gatherv(out, 4, raw(1), integer(1), 0, .comm)
        }

	    else {
            .npRmpi_transport_trace(
                role = "worker",
                event = "worker.exec.send.robj",
                fields = list(tag = tag)
            )
	        mpi.send.Robj(out,0,tag,.comm)
	    }       
	    }
    else {
	        out <- tryCatch(.npRmpi_eval_scmd(scmd.arg$scmd, scmd.arg$arg, envir = parent.frame()),
	                        error = function(e) e)
	    }
    .npRmpi_transport_trace(
        role = "worker",
        event = "worker.exec.done",
        fields = list(tag = tag, ret = ret)
    )
}

mpi.close.Rslaves <- function(dellog=TRUE, comm=1, force=FALSE){
    if (mpi.comm.size(comm) < 2){
    err <-paste("It seems no slaves running on comm", comm)
    stop(err)
    }
    if (isTRUE(getOption("npRmpi.reuse.slaves", FALSE)) && !isTRUE(force)) {
        # Soft-close: keep the slave daemons alive for reuse in this session.
        # This avoids repeated spawn/merge/teardown cycles that can hang/crash
        # on some MPI stacks (notably MPICH on macOS).
        reset.ok <- tryCatch({
          fn <- get(".npRmpi_session_reset_spmd_state_allranks",
                    envir = asNamespace("npRmpi"),
                    mode = "function",
                    inherits = FALSE)
          isTRUE(fn(comm = comm, strict = FALSE, where = "mpi.close.Rslaves() soft-close reset"))
        }, error = function(e) FALSE)
        if (isTRUE(reset.ok))
            return(invisible(0L))
        .np_warning(
            "soft-close SPMD reset failed; falling back to hard-close to avoid stale sequence state",
            call. = FALSE
        )
    }
    # Tell slaves to exit their daemon loop cleanly.
    # The spawned slave daemon (`inst/slavedaemon.R`) treats the character
    # token "kaerb" as a shutdown signal.
    recv.timeout <- tryCatch({
      fn <- get0(".npRmpi_session_recv_timeout", envir = asNamespace("npRmpi"), inherits = FALSE)
      if (is.function(fn)) as.numeric(fn()) else 0
    }, error = function(e) 0)
    if (!is.finite(recv.timeout) || recv.timeout <= 0)
      recv.timeout <- 0

    if (recv.timeout > 0)
      base::setTimeLimit(elapsed = recv.timeout, transient = TRUE)
    shut <- try(mpi.bcast.cmd(cmd="kaerb", rank=0, comm=comm), silent=TRUE)
    if (recv.timeout > 0)
      base::setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
    if (inherits(shut, "try-error")) {
      msg <- as.character(shut)
      .np_warning(
        paste(
          "timed out or failed while broadcasting slave shutdown token.",
          "Proceeding with communicator teardown;",
          "ensure stale mpiexec/slave daemons are cleaned if this recurs.",
          paste(msg, collapse = " ")
        ),
        call. = FALSE
      )
    }
    if (.Platform$OS!="windows"){
        if (dellog && mpi.comm.size(0) < mpi.comm.size(comm)){
        tmp <- paste(Sys.getpid(),"+",comm,sep="")  
        logfile <- paste("*.",tmp,".*.log", sep="")
        logfiles <- Sys.glob(logfile)
        if (length(logfiles) >= 1)
            unlink(logfiles)
        }
    }
#     mpi.barrier(comm)
    if (comm > 0){
        # `comm` is an intracommunicator created via `MPI_Intercomm_merge()`.
        # Using `MPI_Comm_disconnect()` here has been observed to destabilize
        # subsequent spawn/merge cycles with MPICH; `MPI_Comm_free()` is the
        # correct teardown for an intracommunicator.
        if (is.loaded("mpi_comm_free")) {
            res <- mpi.comm.free(comm)
        } else if (is.loaded("mpi_comm_disconnect")) {
            res <- mpi.comm.disconnect(comm)
        } else {
            res <- mpi.comm.free(comm)
        }
    }
#   mpi.comm.set.errhandler(0)
    invisible(res)
}

tailslave.log <- function(nlines=3,comm=1){
    if (mpi.comm.size(comm)==0)
    stop ("It seems no slaves running")
    tmp <- paste(Sys.getpid(),"+",comm,sep="")  
    logfile <- paste("*.",tmp,".*.log", sep="")
    logfiles <- Sys.glob(logfile)
    if (length(logfiles) == 0)
    stop("It seems no slave log files.")
    system(paste("tail -",nlines," ", logfile,sep=""))
}

mpi.apply <- function(X, FUN, ...,  comm=1){
    n <- length(X)
    nslaves <- mpi.comm.size(comm)-1
     if (nslaves < n)
        stop("data length must be at most total slave size")
    if (!is.function(FUN))
        stop("FUN is not a function")
    length(list(...)) #test for any non existing R objects
    tag <- floor(runif(1,1,1000))    
    mpi.bcast.cmd(.mpi.worker.apply, n=n, tag=tag, comm=comm)
    #mpi.bcast(as.integer(c(tag,n)),type=1,comm=comm)
    mpi.bcast.Robj(list(FUN=FUN,dot.arg=list(...)),rank=0,comm=comm)
    if (n < nslaves)
        X=c(X,as.list(integer( nslaves-n)))
    mpi.scatter.Robj(c(list("master"),as.list(X)),root=0,comm=comm)

    out <- as.list(integer(n))
    for (i in seq_len(n)){
       tmp<- mpi.recv.Robj(mpi.any.source(),tag,comm)
       src <- mpi.get.sourcetag()[1]
       out[[src]]<- tmp
    }
    out
}

.mpi.worker.apply <- function(n, tag){
    #assign(".mpi.err", FALSE,  envir = .GlobalEnv)
	.comm <- 1
    #tag.n <- mpi.bcast(integer(2), type=1, comm=.comm)
    #tag <- tag.n[1]
    #n <- tag.n[2]
    tmpfunarg <- mpi.bcast.Robj(rank=0, comm=.comm)
    .tmpfun <- tmpfunarg$FUN
    dotarg <- tmpfunarg$dot.arg
    tmpdata.arg <- list(mpi.scatter.Robj(root=0,comm=.comm))
    if (mpi.comm.rank(.comm) <= n){
        out <- tryCatch(do.call(.tmpfun, c(tmpdata.arg, dotarg)),
                        error = function(e)
                          structure(conditionMessage(e),
                                    class = "try-error",
                                    condition = e))
        mpi.send.Robj(out,0,tag,.comm)
    }
}

mpi.iapply <- function(X, FUN, ...,  comm=1, sleep=0.01){
    n <- length(X)
    nslaves <- mpi.comm.size(comm)-1
     if (nslaves < n)
        stop("data length must be at most total slave size")
    if (!is.function(FUN))
        stop("FUN is not a function")
    length(list(...)) #test for any non existing R objects
    tag <- floor(runif(1,1,1000))    
    mpi.bcast.cmd(.mpi.worker.apply, n=n, tag=tag,comm=comm)
    #mpi.bcast(as.integer(c(tag,n)),type=1,comm=comm)
    mpi.bcast.Robj(list(FUN=FUN,dot.arg=list(...)),rank=0,comm=comm)
    if (n < nslaves)
        X=c(X,as.list(integer( nslaves-n)))
    mpi.scatter.Robj(c(list("master"),as.list(X)),root=0,comm=comm)

    out <- as.list(integer(n))
    done=0
	anysource=mpi.any.source()
    repeat {
       if (mpi.iprobe(anysource,tag,comm)){ 
       srctag <- mpi.get.sourcetag()
       charlen <- mpi.get.count(type=4)
           tmp <- unserialize(mpi.recv(x = raw(charlen), type = 4, srctag[1], 
            srctag[2], comm))
           out[[srctag[1]]]<- tmp
       done=done+1
       }
       if (done < n)
       Sys.sleep(sleep)
       else break
    }
	gc()
    out
}

mpi.parSim <- function(n=100,rand.gen=rnorm, rand.arg=NULL, 
            statistic, nsim=100, run=1, slaveinfo=FALSE, sim.seq=NULL,
            simplify=TRUE, comm=1, ...){
	sim.seq=NULL
    if (mpi.comm.size(comm) < 2)
        stop("It seems no slaves running.")
    if (!is.function(rand.gen))
        stop("rand.gen is not a function")
    if (!is.function(statistic))
        stop("statistic is not a function")
    if (!is.null(rand.arg))
        if (!is.list(rand.arg))
            stop("rand.arg is not a list")
    if (length(list(...))>0)
        deparse(list(...))
        
    slave.num <- mpi.comm.size(comm)-1
    if (!is.null(sim.seq))
        if (!is.integer(sim.seq))
            stop("sim.seq is not an integer vector")
        else if (min(sim.seq)<1 && max(sim.seq)>slave.num && 
                length(sim.seq)!=slave.num*run)
            stop("sim.seq is not in right order")

    mpi.bcast.cmd(.mpi.worker.sim, n=n, nsim=nsim, run=run, comm=comm)  
    mpi.bcast.Robj(list(rand.gen=rand.gen, rand.arg=rand.arg,
                        stat=statistic, stat.arg=list(...)), comm=comm)

    #nnr <- c(n,nsim,run)
    #mpi.bcast(as.integer(nnr),type=1, comm=comm)
    result <- as.list(integer(slave.num*run))

    if (!is.null(sim.seq)){
        for ( i in seq_len(slave.num * run)){
            result[[i]] <- mpi.recv.Robj(source=sim.seq[i], tag=8, comm=comm)
            mpi.send(as.integer(i), type=1, dest=sim.seq[i], tag=88, comm=comm)
        }
        return(.simplify(slave.num*run, result, simplify, nsim))
    }

    i <- 0
    anysrc <- mpi.any.source()
    anytag <- mpi.any.tag()
    mpi.parSim.tmp <- integer(slave.num*run)
    while (i < slave.num*run){
        i <- i+1
        result[[i]] <- mpi.recv.Robj(source=anysrc, tag=8, comm=comm)
        src <- mpi.get.sourcetag()[1]
        mpi.send(as.integer(i), type=1, dest=src, tag=88, comm=comm)
        mpi.parSim.tmp[i] <- src
    }
    if (slaveinfo){
        slavename <- paste("slave", seq_len(slave.num), sep="")
        cat("Finished slave jobs summary:\n")
        for (i in seq_len(slave.num)){
            if (i < 10)
                cat(slavename[i], " finished",sum(mpi.parSim==i), "job(s)\n")
            else
                cat(slavename[i], "finished",sum(mpi.parSim==i), "job(s)\n")
        }
    }
	#assign(".mpi.parSim", mpi.parSim.tmp,  envir = .GlobalEnv)
    .simplify(slave.num*run, result, simplify, nsim)
}

.mpi.worker.sim <- function(n, nsim, run){
	.comm <- 1
    tmpdata <- mpi.bcast.Robj(comm=.comm)
    rand.arg <- tmpdata$rand.arg
    stat.arg <- tmpdata$stat.arg
    
    .tmp.rand.gen <- tmpdata$rand.gen
    .tmp.statistic <- tmpdata$stat
    
    #nnr <- mpi.bcast(integer(3), type=1, comm=.comm)
    #n <- nnr[1];  nsim <- nnr[2];  run <- nnr[3]

    i <- 0
    slave.num <- mpi.comm.size(.comm)-1
    
    while( i < slave.num*(run-1)+1){
        out <- replicate(nsim, do.call(.tmp.statistic, c(list(do.call(.tmp.rand.gen, 
                            c(list(n),rand.arg))), stat.arg)))
       
        mpi.send.Robj(obj=out, dest=0, tag=8, comm=.comm)
        i <- mpi.recv(integer(1), type=1, source=0, tag=88, comm=.comm)
    }
}

#from snow
.docall <- function(fun, args, envir = parent.frame()) {
    if ((is.character(fun) && length(fun) == 1) || is.name(fun)) {
        fname <- as.character(fun)
        fun <- get0(fname, envir = envir, mode = "function", inherits = TRUE)
        if (is.null(fun))
            fun <- get0(fname, envir = .GlobalEnv, mode = "function", inherits = FALSE)
        if (is.null(fun))
            stop(sprintf("object '%s' of mode 'function' was not found", fname), call. = FALSE)
    }
    enquote <- function(x) as.call(list(as.name("quote"), x))
    do.call(fun, lapply(args, enquote))
}

.splitIndices <- function(nx, ncl) {
    #i <- 1:nx;
    #structure(split(i, cut(i, ncl)), names=NULL)
    x <- seq_len(nx)
    r <- nx/ncl
    ii <- seq.int(0L, ncl - 1L) * r
    if (nx < ncl)
        intv <- seq.int(0L, ncl)
    else
        intv <- c(x[round(1 + ii)] - 1L, nx)
    structure(split(x, cut(x, intv)), names = NULL)
}

mpi.parMM <- function(A, B, job.num=mpi.comm.size(comm)-1, comm=1){
    splitRows <- function(x, ncl)
        lapply(.splitIndices(nrow(x), ncl), function(i) x[i,, drop=FALSE])    
    .docall(rbind, mpi.applyLB(splitRows(A, job.num), 
        get("%*%"), B, comm=comm))
}
    
mpi.iparMM <- function(A, B, comm=1, sleep=0.01){
    splitRows <- function(x, ncl)
        lapply(.splitIndices(nrow(x), ncl), function(i) x[i,, drop=FALSE])    
    .docall(rbind, mpi.iapply(splitRows(A, mpi.comm.size(comm)-1), 
        get("%*%"), B, comm=comm, sleep=sleep))
}    

mpi.applyLB <- function(X, FUN, ...,  apply.seq=NULL, comm=1){
	apply.seq=NULL
    n <- length(X)
    slave.num <- mpi.comm.size(comm)-1
    .npRmpi_transport_trace(
        role = "master",
        event = "applylb.start",
        fields = list(n = n, slave_num = slave.num, comm = comm)
    )
    if (slave.num < 1)
        stop("There are no slaves running")
    if (n <= slave.num) {
        .npRmpi_clear_applylb_cache()
        return (mpi.apply(X,FUN,...,comm=comm))
    }    
    if (!is.function(FUN))
        stop("FUN is not a function")
    length(list(...))
    if (!is.null(apply.seq))
        if (!is.integer(apply.seq))
            stop("apply.seq is not an integer vector")
        else if (min(apply.seq)<1 && max(apply.seq)>slave.num && 
                length(apply.seq)!=n)
            stop("apply.seq is not in right order")
            
    mpi.bcast.cmd(.mpi.worker.applyLB, n=n, comm=comm)
    #mpi.bcast(as.integer(n),type=1,comm=comm)
    mpi.bcast.Robj(list(FUN=FUN,dot.arg=list(...)),rank=0,comm=comm)
    out <- as.list(integer(n))
    mpi.anysource <- mpi.any.source()
    mpi.anytag <- mpi.any.tag()
    for (i in seq_len(slave.num))
        mpi.send.Robj(list(data.arg=list(X[[i]])), dest=i,tag=i, comm=comm)
  
    if (!is.null(apply.seq)){
        for ( i in seq_len(n)){
            tmp <- mpi.recv.Robj(source=apply.seq[i], tag=mpi.anytag, comm=comm)
            tag <- mpi.get.sourcetag()[2]
            out[[tag]]<- tmp
            j <- i+slave.num
            if (j <= n)
                mpi.send.Robj(list(data.arg=list(X[[j]])), dest=apply.seq[i],tag=j, comm=comm)
            else
                mpi.send.Robj(as.integer(0),dest=apply.seq[i],tag=j,comm=comm)
        }
        .npRmpi_transport_trace(
            role = "master",
            event = "applylb.done",
            fields = list(n = n, slave_num = slave.num, comm = comm, apply_seq = TRUE)
        )
        return(out)
    }
    # .mpi.applyLB <- integer(n)
	mpi.seq.tmp <- integer(n)
    for (i in seq_len(n)){
       tmp<- mpi.recv.Robj(mpi.anysource,mpi.anytag,comm)
       srctag <- mpi.get.sourcetag()
       out[[srctag[2]]]<- tmp
       mpi.seq.tmp[i] <- srctag[1]
       j <- i+slave.num
       if (j <= n)
            mpi.send.Robj(list(data.arg=list(X[[j]])), dest=srctag[1],tag=j, comm=comm)
       else
            mpi.send.Robj(as.integer(0),dest=srctag[1],tag=j,comm=comm)
    }
 	#assign(".mpi.applyLB",mpi.seq.tmp, envir = .GlobalEnv)
    .npRmpi_transport_trace(
        role = "master",
        event = "applylb.done",
        fields = list(n = n, slave_num = slave.num, comm = comm, apply_seq = FALSE)
    )
	out
}

.mpi.worker.applyLB <- function(n){
    #assign(".mpi.err", FALSE,  envir = .GlobalEnv)
	.comm <- 1
    .npRmpi_transport_trace(
        role = "worker",
        event = "worker.applylb.start",
        fields = list(n = n, comm = .comm)
    )
    #n <- mpi.bcast(integer(1), type=1, comm=.comm)
    tmpfunarg <- mpi.bcast.Robj(rank=0, comm=.comm)
    .tmpfun <- tmpfunarg$FUN
    dotarg <- tmpfunarg$dot.arg
    mpi.anytag <- mpi.any.tag()
    repeat {
        tmpmsg <- mpi.recv.Robj(source=0,tag=mpi.anytag, comm=.comm)
        tag <- mpi.get.sourcetag()[2]
        if (tag > n)
            break
        tmpdata.arg <- if (is.list(tmpmsg)) tmpmsg$data.arg else NULL
        if (is.null(tmpdata.arg)) {
            .npRmpi_transport_trace(
                role = "worker",
                event = "worker.applylb.malformed_task",
                fields = list(tag = tag)
            )
            out <- structure(
                "mpi.applyLB worker received malformed task payload",
                class = "try-error"
            )
            mpi.send.Robj(out,0,tag,.comm)
            next
        }
        out <- tryCatch(do.call(.tmpfun, c(tmpdata.arg, dotarg)),
                        error = function(e)
                          structure(conditionMessage(e),
                                    class = "try-error",
                                    condition = e))
        #if (.mpi.err)
        #    print(geterrmessage())
        mpi.send.Robj(out,0,tag,.comm)
    }
    .npRmpi_transport_trace(
        role = "worker",
        event = "worker.applylb.done",
        fields = list(n = n, comm = .comm)
    )
}

mpi.iapplyLB <- function(X, FUN, ...,  apply.seq=NULL, comm=1, sleep=0.01){
	apply.seq=NULL
    n <- length(X)
    slave.num <- mpi.comm.size(comm)-1
    if (slave.num < 1)
        stop("There are no slaves running")
    if (n <= slave.num) {
        .npRmpi_clear_applylb_cache()
        return (mpi.iapply(X,FUN,...,comm=comm,sleep=sleep))
    }
    if (!is.function(FUN))
        stop("FUN is not a function")
    if (slave.num > 2000)
        stop("Total slaves are more than nonblock send/receive can handle")
    length(list(...))
    if (!is.null(apply.seq))
        if (!is.integer(apply.seq))
            stop("apply.seq is not an integer vector")
        else if (min(apply.seq)<1 && max(apply.seq)>slave.num &&
                length(apply.seq)!=n)
            stop("apply.seq is not in right order")

    mpi.bcast.cmd(.mpi.worker.applyLB, n=n, comm=comm)
    #mpi.bcast(as.integer(n),type=1,comm=comm)
    mpi.bcast.Robj(list(FUN=FUN,dot.arg=list(...)),rank=0,comm=comm)
    out <- as.list(integer(n))
    mpi.anysource <- mpi.any.source()
    mpi.anytag <- mpi.any.tag()
    for (i in seq_len(slave.num))
        mpi.send.Robj(list(data.arg=list(X[[i]])), dest=i,tag=i,comm=comm)
    #for (i in 1:slave.num)
    #    mpi.waitany(slave.num)

    if (!is.null(apply.seq)){
       i=0
       repeat {
        if (mpi.iprobe(apply.seq[i+1],mpi.anytag,comm)){
            i=i+1
            j <- i+slave.num
            if ( j <= n)
                mpi.send.Robj(list(data.arg=list(X[[j]])), dest=apply.seq[i],tag=j, comm=comm) 
            else
                mpi.send.Robj(as.integer(0),dest=apply.seq[i],tag=j,comm=comm)  
            charlen <- mpi.get.count(type=4)
            tag <- mpi.get.sourcetag()[2]
            tmp <- unserialize(mpi.recv(x = raw(charlen), type = 4, apply.seq[i], tag, comm))
            out[[tag]]<- tmp
            #mpi.wait(0)
        }
      if (i < n)
         Sys.sleep(sleep)
      else break
      }
      return(out)
    }
    mpi.seq.tmp <- integer(n)
    i=0
    repeat {
        if (mpi.iprobe(mpi.anysource,mpi.anytag,comm)){
            i=i+1
            srctag <- mpi.get.sourcetag()
            src <- srctag[1]
            tag <- srctag[2]
            j <- i+slave.num
            if ( j <= n)
                mpi.send.Robj(list(data.arg=list(X[[j]])), dest=src,tag=j, comm=comm)
            else
                mpi.send.Robj(as.integer(0),dest=src,tag=j,comm=comm)
            charlen <- mpi.get.count(type=4)
            tmp <- unserialize(mpi.recv(x = raw(charlen), type = 4, src, tag, comm))
            out[[tag]]<- tmp
            mpi.seq.tmp[i] <- src
            #mpi.wait(src-1)
        }
        if (i < n)
            Sys.sleep(sleep)
        else
            break
    }
  	#assign(".mpi.applyLB",mpi.seq.tmp, envir = .GlobalEnv)
	gc()
    out
}

#.mpi.worker.iapplyLB <- function(){
#    assign(".mpi.err", FALSE,  envir = .GlobalEnv)
#    n <- mpi.bcast(integer(1), type=1, comm=.comm)
#    tmpfunarg <- mpi.bcast.Robj(rank=0, comm=.comm)
#    .tmpfun <- tmpfunarg$fun
#    dotarg <- tmpfunarg$dot.arg
#    mpi.anytag <- mpi.any.tag()
#    repeat {
#        tmpdata.arg <- mpi.recv.Robj(source=0,tag=mpi.anytag, comm=.comm)$data.arg
#        tag <- mpi.get.sourcetag()[2]
#        if (tag > n)
#            break
#        out <- try(do.call(.tmpfun, c(tmpdata.arg, dotarg)),TRUE)
#        #if (.mpi.err)
#        #    print(geterrmessage())
#        mpi.wait(0)
#        mpi.isend.Robj(out,0,tag,.comm)
#    }
#    mpi.wait(0)
#}

.simplify <- function(n, answer, simplify, len=1, recursive=FALSE){
    if (simplify && length(answer)&&length(common.len <- unique(unlist(lapply(answer, 
        length)))) == 1 ) {
        if (common.len == len) 
            unlist(answer, recursive = recursive)
        else if (common.len > len) 
            array(unlist(answer, recursive = recursive),
                  dim = c(common.len/len, n*len),
                  dimnames = list(names(answer[[1]]), names(answer)))
        else answer
    }
    else answer
}

mpi.parLapply <- function(X, FUN, ..., job.num=mpi.comm.size(comm)-1, apply.seq=NULL, comm=1){
    if (job.num < 2)
        stop("job.num is at least 2.")
    splitList <- function(x, ncl)
        lapply(.splitIndices(length(X), ncl), function(i) X[i])
    .docall(c, mpi.applyLB(splitList(X, job.num), 
        lapply, FUN, ..., apply.seq=apply.seq, comm=comm))
}

mpi.iparLapply <- function(X, FUN, ..., job.num=mpi.comm.size(comm)-1, apply.seq=NULL, comm=1, sleep=0.01){
    if (job.num < 2)
        stop("job.num is at least 2.")
    splitList <- function(X, ncl)
        lapply(.splitIndices(length(X), ncl), function(i) X[i])
    .docall(c, mpi.iapplyLB(splitList(X, job.num), 
        lapply, FUN, ..., apply.seq=apply.seq, comm=comm, sleep=sleep))
}

mpi.parSapply <- function (X, FUN, ..., job.num=mpi.comm.size(comm)-1, apply.seq=NULL, 
                simplify = TRUE, USE.NAMES = TRUE, comm=1) 
{
    FUN <- match.fun(FUN)
    answer <- mpi.parLapply(as.list(X),FUN,...,job.num=job.num,apply.seq=apply.seq,comm=comm)
    if (USE.NAMES && is.character(X) && is.null(names(answer))) 
        names(answer) <- X
    .simplify(length(X),answer, simplify)
}

mpi.iparSapply <- function (X, FUN, ..., job.num=mpi.comm.size(comm)-1, apply.seq=NULL, 
                simplify = TRUE, USE.NAMES = TRUE, comm=1,sleep=0.01) 
{
    FUN <- match.fun(FUN)
    answer <- mpi.iparLapply(as.list(X),FUN,...,job.num=job.num,apply.seq=apply.seq,comm=comm,sleep=sleep)
    if (USE.NAMES && is.character(X) && is.null(names(answer))) 
        names(answer) <- X
    .simplify(length(X),answer, simplify)
}

.mpi_make_replicate_fun <- function(expr, env = parent.frame()) {
    as.function(c(alist(... = ), expr), envir = env)
}

mpi.parReplicate <- function(n,  expr, job.num=mpi.comm.size(comm)-1, apply.seq=NULL,
                                simplify = TRUE, comm=1){
    expr_fun <- .mpi_make_replicate_fun(substitute(expr), env = parent.frame())
    mpi.parSapply(integer(n), expr_fun, 
    job.num=job.num, apply.seq=apply.seq, simplify = simplify, comm=comm)
}

mpi.iparReplicate <- function(n,  expr, job.num=mpi.comm.size(comm)-1, apply.seq=NULL,
                                simplify = TRUE, comm=1,sleep=0.01){
    expr_fun <- .mpi_make_replicate_fun(substitute(expr), env = parent.frame())
    mpi.iparSapply(integer(n), expr_fun, 
    job.num=job.num, apply.seq=apply.seq, simplify = simplify, comm=comm,sleep=sleep)
}

mpi.parRapply <- function(X,FUN,...,job.num=mpi.comm.size(comm)-1,apply.seq=NULL,comm=1){
    if (job.num < 2)
        stop("job.num is at least 2.")
    splitRows <- function(X, ncl)
        lapply(.splitIndices(nrow(X), ncl), function(i) X[i,, drop=FALSE])
    .docall(c, mpi.applyLB(splitRows(X,job.num), apply, 1, FUN, ..., 
            apply.seq=apply.seq, comm=comm))
}

mpi.iparRapply <- function(X,FUN,...,job.num=mpi.comm.size(comm)-1,apply.seq=NULL,comm=1,sleep=0.01){
    if (job.num < 2)
        stop("job.num is at least 2.")
    splitRows <- function(X, ncl)
        lapply(.splitIndices(nrow(X), ncl), function(i) X[i,, drop=FALSE])
    .docall(c, mpi.iapplyLB(splitRows(X,job.num), apply, 1, FUN, ..., 
            apply.seq=apply.seq, comm=comm,sleep=sleep))
}

mpi.parCapply <- function(X,FUN,...,job.num=mpi.comm.size(comm)-1,apply.seq=NULL,comm=1){
    if (job.num < 2)
        stop("job.num is at least 2.")
    splitCols <- function(X, ncl)
        lapply(.splitIndices(ncol(X), ncl), function(i) X[,i, drop=FALSE])
    .docall(c, mpi.applyLB(splitCols(X,job.num), apply, 2, FUN, ..., 
            apply.seq=apply.seq, comm=comm))
}

mpi.iparCapply <- function(X,FUN,...,job.num=mpi.comm.size(comm)-1,apply.seq=NULL,comm=1,sleep=0.01){
    if (job.num < 2)
        stop("job.num is at least 2.")
    splitCols <- function(X, ncl)
        lapply(.splitIndices(ncol(X), ncl), function(i) X[,i, drop=FALSE])
    .docall(c, mpi.iapplyLB(splitCols(X,job.num), apply, 2, FUN, ..., 
            apply.seq=apply.seq, comm=comm,sleep=sleep))
}

mpi.parApply <- function(X, MARGIN, FUN, ..., job.num = mpi.comm.size(comm)-1,
                    apply.seq=NULL, comm=1)
{
    FUN <- match.fun(FUN)
    d <- dim(X)
    dl <- length(d)
    if(dl == 0)
    stop("dim(X) must have a positive length")
    ds <- seq_len(dl)

    if(length(oldClass(X)) > 0)
    X <- if(dl == 2) as.matrix(X) else as.array(X)
    dn <- dimnames(X)

    s.call <- ds[-MARGIN]
    s.ans  <- ds[MARGIN]
    d.call <- d[-MARGIN]
    d.ans  <- d[MARGIN]
    dn.call<- dn[-MARGIN]
    dn.ans <- dn[MARGIN]

    d2 <- prod(d.ans)
    if(d2 == 0) {
        newX <- array(vector(typeof(X), 1), dim = c(prod(d.call), 1))
        ans <- FUN(if(length(d.call) < 2) newX[,1] else
                   array(newX[,1], d.call, dn.call), ...)
        return(if(is.null(ans)) ans else if(length(d.call) < 2) ans[1][-1]
               else array(ans, d.ans, dn.ans))
    }
    newX <- aperm(X, c(s.call, s.ans))
    dim(newX) <- c(prod(d.call), d2)
    if(length(d.call) < 2) {
        if (length(dn.call)) dimnames(newX) <- c(dn.call, list(NULL))
        ans <- mpi.parLapply(seq_len(d2), function(i, ...) FUN(newX[,i], ...), 
                ..., job.num = job.num, apply.seq=apply.seq, comm=comm )
    } 
    else
        ans <- mpi.parLapply(seq_len(d2),
            function(i, ...) FUN(array(newX[,i], d.call, dn.call), ...),
                ..., job.num = job.num, apply.seq=apply.seq, comm=comm)
                    
    ans.list <- is.recursive(ans[[1]])
    l.ans <- length(ans[[1]])

    ans.names <- names(ans[[1]])
    if(!ans.list)
        ans.list <- any(unlist(lapply(ans, length)) != l.ans)
    if(!ans.list && length(ans.names)) {
        all.same <- sapply(ans, function(X) identical(names(X), ans.names))
        if (!all(all.same)) ans.names <- NULL
    }
    len.a <- if(ans.list) d2 else length(ans <- unlist(ans, recursive = FALSE))
    if(length(MARGIN) == 1 && len.a == d2) {
        names(ans) <- if(length(dn.ans[[1]])) dn.ans[[1]]
        return(ans)
    }
    if(len.a == d2)
        return(array(ans, d.ans, dn.ans))
    if(len.a > 0 && len.a %% d2 == 0)
        return(array(ans, c(len.a %/% d2, d.ans),
                     if(is.null(dn.ans)) {
                         if(!is.null(ans.names)) list(ans.names,NULL)
                     } else c(list(ans.names), dn.ans)))
    return(ans)
}

mpi.iparApply <- function(X, MARGIN, FUN, ..., job.num = mpi.comm.size(comm)-1,
                    apply.seq=NULL, comm=1,sleep=0.01)
{
    FUN <- match.fun(FUN)
    d <- dim(X)
    dl <- length(d)
    if(dl == 0)
    stop("dim(X) must have a positive length")
    ds <- seq_len(dl)

    if(length(oldClass(X)) > 0)
    X <- if(dl == 2) as.matrix(X) else as.array(X)
    dn <- dimnames(X)

    s.call <- ds[-MARGIN]
    s.ans  <- ds[MARGIN]
    d.call <- d[-MARGIN]
    d.ans  <- d[MARGIN]
    dn.call<- dn[-MARGIN]
    dn.ans <- dn[MARGIN]

    d2 <- prod(d.ans)
    if(d2 == 0) {
        newX <- array(vector(typeof(X), 1), dim = c(prod(d.call), 1))
        ans <- FUN(if(length(d.call) < 2) newX[,1] else
                   array(newX[,1], d.call, dn.call), ...)
        return(if(is.null(ans)) ans else if(length(d.call) < 2) ans[1][-1]
               else array(ans, d.ans, dn.ans))
    }
    newX <- aperm(X, c(s.call, s.ans))
    dim(newX) <- c(prod(d.call), d2)
    if(length(d.call) < 2) {
        if (length(dn.call)) dimnames(newX) <- c(dn.call, list(NULL))
        ans <- mpi.iparLapply(seq_len(d2), function(i, ...) FUN(newX[,i], ...), 
                ..., job.num = job.num, apply.seq=apply.seq, comm=comm, sleep=sleep )
    } 
    else
        ans <- mpi.iparLapply(seq_len(d2),
            function(i, ...) FUN(array(newX[,i], d.call, dn.call), ...),
                ..., job.num = job.num, apply.seq=apply.seq, comm=comm, sleep=sleep)
                    
    ans.list <- is.recursive(ans[[1]])
    l.ans <- length(ans[[1]])

    ans.names <- names(ans[[1]])
    if(!ans.list)
        ans.list <- any(unlist(lapply(ans, length)) != l.ans)
    if(!ans.list && length(ans.names)) {
        all.same <- sapply(ans, function(X) identical(names(X), ans.names))
        if (!all(all.same)) ans.names <- NULL
    }
    len.a <- if(ans.list) d2 else length(ans <- unlist(ans, recursive = FALSE))
    if(length(MARGIN) == 1 && len.a == d2) {
        names(ans) <- if(length(dn.ans[[1]])) dn.ans[[1]]
        return(ans)
    }
    if(len.a == d2)
        return(array(ans, d.ans, dn.ans))
    if(len.a > 0 && len.a %% d2 == 0)
        return(array(ans, c(len.a %/% d2, d.ans),
                     if(is.null(dn.ans)) {
                         if(!is.null(ans.names)) list(ans.names,NULL)
                     } else c(list(ans.names), dn.ans)))
    return(ans)
}

#mpi.parallel.sim <- mpi.parSim
