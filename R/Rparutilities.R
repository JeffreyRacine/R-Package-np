### Copyright (C) 2002 Hao Yu 

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
    if (!mpi.is.master())
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
        ranks <- 1:(size-1)
        commm <- paste(comm, ")",sep="")
        if (size > 10){
        rank0 <- paste("master  (rank 0 , comm", commm)
            ranks <- c(paste(ranks[1:9]," ",sep=""), ranks[10:(size-1)])
        }
        else
        rank0 <- paste("master (rank 0, comm", commm)
        cat(rank0, "of size", size, "is running on:",master, "\n")
        slavename <- paste("slave", ranks,sep="")
        ranks <- paste("(rank ",ranks, ", comm ",slavecomm,")", sep="")
		if (short && size > 8){
          for (i in 1:3) {
            cat(slavename[i], ranks[i], "of size",size, 
          "is running on:",slavehost[i], "\n")	
		  }
		  cat("... ... ...\n")
		  for (i in (size-2):(size-1)){
		    cat(slavename[i], ranks[i], "of size",size, 
          "is running on:",slavehost[i], "\n")
		  }
		}
		else {
          for (i in 1:(size-1)){
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
    nn <- 0:(length(hosts)-1)
        names(nn) <- base
    nn
}

mpichhosts <- function(){
    if (.Platform$OS != "windows") 
        stop("mpichhosts runs only under MPICH2 for Windows")
    hosts <- system("smpd -get hosts", intern=TRUE)
    
    if (length(hosts) == 0) 
	  hostnames <- "localhost"
	else {
		if (hosts=="default")
		    hostnames <- "localhost"
    	else	{ 
        	hosts <- unlist(strsplit(hosts, " "))
        	hosts <- hosts[which(hosts != "")]
        	hostnames <- NULL
        	for (host in hosts) {
            	hostsmp <- unlist(strsplit(host, ":"))
            	smp <- ifelse(is.na(hostsmp[2]), 1, as.integer(hostsmp[2]))
            	hostbase <- unlist(strsplit(hostsmp[1], "\\."))[1]
            	hostnames <- c(hostnames, rep(hostbase, smp))
        	}
    	}
	}
    base <- "master"
    if (length(hostnames) == 1) {
	  out=0
	  repeat {
        	out <- out + 1
            cpus = length(.Call("RegQuery", as.integer(3), 
            paste("HARDWARE\\DESCRIPTION\\System\\CentralProcessor", 
                   out, sep = "\\"), PACKAGE = "npRmpi"))
            if (cpus == 0) 
                break
        }
	  hostnames=c(hostnames,rep("localhost",out))
    }
    else
	  hostnames = c("localhost", hostnames)
    base <- c(base, paste("slave", 1:(length(hostnames)-1), sep = ""))
    names(hostnames) <- base
    hostnames
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
    if (mpi.comm.size(comm) > 0){
         err <-paste("It seems there are some slaves running on comm ", comm)
         stop(err)
    }
    if (.Platform$OS=="windows"){
        workdrive <- unlist(strsplit(getwd(),":"))[1]
		workdir <- unlist(strsplit(getwd(),"/"))
		if (length(workdir) > 1)
			workdir <-paste(workdir, collapse="\\")
		else
        	workdir <- paste(workdir,"\\")
        localhost <- Sys.getenv("COMPUTERNAME")
        networkdrive <-.Call("RegQuery", as.integer(2),paste("NETWORK\\",workdrive,sep=""), 
                        PACKAGE = "npRmpi")
        remotepath <-networkdrive[which(networkdrive=="RemotePath")+1]
        mapdrive <- as.logical(mapdrive && !is.null(remotepath))
        arg <- c(Rscript, R.home(), workdrive, workdir, localhost, mapdrive, remotepath)
		if (.Platform$r_arch == "i386") 
            realns <- mpi.comm.spawn(slave = system.file("Rslaves32.cmd", 
				package = "npRmpi"), slavearg = arg, nslaves = nslaves, 
				info = 0, root = root, intercomm = intercomm, quiet = quiet)
		else 
			realns <- mpi.comm.spawn(slave = system.file("Rslaves64.cmd", 
			            package = "npRmpi"), slavearg = arg, nslaves = nslaves, 
				info = 0, root = root, intercomm = intercomm, quiet = quiet)
    }
    else{
        tmp <- paste(Sys.getpid(), "+", comm, sep="")   
        if (needlog)
            arg <- c(Rscript, tmp, "needlog", R.home())
        else
            arg <- c(Rscript, tmp , "nolog", R.home())  
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
        if (!quiet) slave.hostinfo(comm)
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
    #if (length(arg) > 0) 
    #    deparse(arg)
    #tag.ret <- c(tag, ret, simplify)
    mpi.bcast.cmd(.mpi.worker.exec, tag=tag, ret=ret, simplify=simplify, comm = comm)
    #mpi.bcast(as.integer(tag.ret), type=1, comm=comm)
    mpi.bcast.Robj(list(scmd=scmd, arg=arg), comm=comm)

    if (ret){
        size <- mpi.comm.size(comm) 
        allcode <- mpi.allgather(integer(2), 1, integer(2*size), comm)
    type <- allcode[seq(3,2*size,2)]
    len <- allcode[seq(4,2*size,2)]
    eqlen <- all(len==len[1])
    if (all(type==1)){
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
                names(out) <- paste("slave",1:(size-1), sep="")
                for (i in 1:(size-1))
            out[[i]]<- out1[lowlen[i]:uplen[i]]
        }
    }
    else if (all(type==2)){
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
                names(out) <- paste("slave",1:(size-1), sep="")
                for (i in 1:(size-1))
            out[[i]]<- out1[lowlen[i]:uplen[i]]
        }
    }
 else if (all(type==4)){
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
                names(out) <- paste("slave",1:(size-1), sep="")
                for (i in 1:(size-1))
                    out[[i]]<- out1[lowlen[i]:uplen[i]]
        }
    }

    else {
            out <- as.list(integer(size-1))
            names(out) <- paste("slave",1:(size-1), sep="")
            for (i in 1:(size-1)){
        tmp<- mpi.recv.Robj(mpi.any.source(),tag,comm)
        src <- mpi.get.sourcetag()[1] 
        out[[src]]<- tmp 
        }
    }
        out
    }
}

.typeindex <- function (x) {
    if(class(x)=="integer")
            as.integer(c(1,length(x)))
    else if (class(x)=="numeric")
            as.integer(c(2,length(x)))
    else if (class(x)=="raw")
            as.integer(c(4,length(x)))

    else
            as.integer(-1)
}

.mpi.worker.exec <- function(tag, ret, simplify){
    #assign(".mpi.err", FALSE,  envir = .GlobalEnv)
    assign(".mpi.err", FALSE)
	.comm <- 1
    #tag.ret <- mpi.bcast(integer(3), type=1, comm=.comm)
    #tag <- tag.ret[1]
    #ret <- as.logical(tag.ret[2])
    #simplify <- as.logical(tag.ret[3])
    scmd.arg <- mpi.bcast.Robj(comm=.comm)

    if (ret){
    size <- mpi.comm.size(.comm)
    myerrcode <- as.integer(0)
    if (length(scmd.arg$arg)>0)
        out <- try(do.call(as.character(scmd.arg$scmd), scmd.arg$arg, envir=.GlobalEnv),TRUE)
    else
        out <- try(eval(scmd.arg$scmd, envir=sys.parent()), TRUE)
    
    if (get(".mpi.err")){
        print(geterrmessage())
        type <- integer(2)
    }
        else {
        type <- .typeindex(out)
        if (is.na(type[2]))
            type[2] <- as.integer(0)    
        }
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
        mpi.send.Robj(out,0,tag,.comm)
    }       
    }
    else {
    if (length(scmd.arg$arg)>0)
            out <- try(do.call(as.character(scmd.arg$scmd), scmd.arg$arg))
        else
            out <- try(eval(scmd.arg$scmd))  
    }
}

mpi.close.Rslaves <- function(dellog=TRUE, comm=1){
    if (mpi.comm.size(comm) < 2){
    err <-paste("It seems no slaves running on comm", comm)
    stop(err)
    }
    mpi.bcast.cmd(cmd=break, rank=0, comm=comm)
    if (.Platform$OS!="windows"){
        if (dellog && mpi.comm.size(0) < mpi.comm.size(comm)){
        tmp <- paste(Sys.getpid(),"+",comm,sep="")  
        logfile <- paste("*.",tmp,".*.log", sep="")
        if (length(system(paste("ls", logfile),TRUE,ignore.stderr=TRUE) )>=1)
            system(paste("rm", logfile))
        }
    }
#     mpi.barrier(comm)
    if (comm >0){
        if (is.loaded("mpi_comm_disconnect"))
            mpi.comm.disconnect(comm) 
        else
            mpi.comm.free(comm)
    }
#   mpi.comm.set.errhandler(0)
}

tailslave.log <- function(nlines=3,comm=1){
    if (mpi.comm.size(comm)==0)
    stop ("It seems no slaves running")
    tmp <- paste(Sys.getpid(),"+",comm,sep="")  
    logfile <- paste("*.",tmp,".*.log", sep="")
    if (length(system(paste("ls", logfile),TRUE,ignore.stderr=TRUE))==0)
    stop("It seems no slave log files.")
    system(paste("tail -",nlines," ", logfile,sep=""))
}

mpi.apply <- function(x, fun, ...,  comm=1){
    n <- length(x)
    nslaves <- mpi.comm.size(comm)-1
     if (nslaves < n)
        stop("data length must be at most total slave size")
    if (!is.function(fun))
        stop("fun is not a function")
    length(list(...)) #test for any non existing R objects
    tag <- floor(runif(1,1,1000))    
    mpi.bcast.cmd(.mpi.worker.apply, n=n, tag=tag, comm=comm)
    #mpi.bcast(as.integer(c(tag,n)),type=1,comm=comm)
    mpi.bcast.Robj(list(fun=fun,dot.arg=list(...)),rank=0,comm=comm)
    if (n < nslaves)
        x=c(x,as.list(integer( nslaves-n)))
    mpi.scatter.Robj(c(list("master"),as.list(x)),root=0,comm=comm)

    out <- as.list(integer(n))
    for (i in 1:n){
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
    .tmpfun <- tmpfunarg$fun
    dotarg <- tmpfunarg$dot.arg
    tmpdata.arg <- list(mpi.scatter.Robj(root=0,comm=.comm))
    if (mpi.comm.rank(.comm) <= n){
        out <- try(do.call(".tmpfun", c(tmpdata.arg, dotarg)),TRUE)
        mpi.send.Robj(out,0,tag,.comm)
    }
}

mpi.iapply <- function(x, fun, ...,  comm=1, sleep=0.01){
    n <- length(x)
    nslaves <- mpi.comm.size(comm)-1
     if (nslaves < n)
        stop("data length must be at most total slave size")
    if (!is.function(fun))
        stop("fun is not a function")
    length(list(...)) #test for any non existing R objects
    tag <- floor(runif(1,1,1000))    
    mpi.bcast.cmd(.mpi.worker.apply, n=n, tag=tag,comm=comm)
    #mpi.bcast(as.integer(c(tag,n)),type=1,comm=comm)
    mpi.bcast.Robj(list(fun=fun,dot.arg=list(...)),rank=0,comm=comm)
    if (n < nslaves)
        x=c(x,as.list(integer( nslaves-n)))
    mpi.scatter.Robj(c(list("master"),as.list(x)),root=0,comm=comm)

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
        for ( i in 1:(slave.num*run)){
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
        slavename <- paste("slave",1:slave.num, sep="")
        cat("Finished slave jobs summary:\n")
        for (i in 1:slave.num){
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
        out <- replicate(nsim, do.call(".tmp.statistic", c(list(do.call(".tmp.rand.gen", 
                            c(list(n),rand.arg))), stat.arg)))
       
        mpi.send.Robj(obj=out, dest=0, tag=8, comm=.comm)
        i <- mpi.recv(integer(1), type=1, source=0, tag=88, comm=.comm)
    }
}

#from snow
.docall <- function(fun, args) {
    if ((is.character(fun) && length(fun) == 1) || is.name(fun))
        fun <- get(as.character(fun), envir = .GlobalEnv, mode = "function")
    enquote <- function(x) as.call(list(as.name("quote"), x))
    do.call("fun", lapply(args, enquote))
}

.splitIndices <- function(nx, ncl) {
    #i <- 1:nx;
    #structure(split(i, cut(i, ncl)), names=NULL)
    x <- 1:nx
    r <- nx/ncl
    ii <- 0:(ncl - 1) * r
    if (nx < ncl)
        intv <- 0:ncl
    else
        intv <- c(x[round(1 + ii)]-1,nx)
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

mpi.applyLB <- function(x, fun, ...,  apply.seq=NULL, comm=1){
	apply.seq=NULL
    n <- length(x)
    slave.num <- mpi.comm.size(comm)-1
    if (slave.num < 1)
        stop("There are no slaves running")
    if (n <= slave.num) {
        if (exists(".mpi.applyLB")) 
			rm(".mpi.applyLB",  envir=.GlobalEnv)
        return (mpi.apply(x,fun,...,comm=comm))
    }    
    if (!is.function(fun))
        stop("fun is not a function")
    length(list(...))
    if (!is.null(apply.seq))
        if (!is.integer(apply.seq))
            stop("apply.seq is not an integer vector")
        else if (min(apply.seq)<1 && max(apply.seq)>slave.num && 
                length(apply.seq)!=n)
            stop("apply.seq is not in right order")
            
    mpi.bcast.cmd(.mpi.worker.applyLB, n=n, comm=comm)
    #mpi.bcast(as.integer(n),type=1,comm=comm)
    mpi.bcast.Robj(list(fun=fun,dot.arg=list(...)),rank=0,comm=comm)
    out <- as.list(integer(n))
    mpi.anysource <- mpi.any.source()
    mpi.anytag <- mpi.any.tag()
    for (i in 1:slave.num)
        mpi.send.Robj(list(data.arg=list(x[[i]])), dest=i,tag=i, comm=comm)
  
    if (!is.null(apply.seq)){
        for ( i in 1:n){
            tmp <- mpi.recv.Robj(source=apply.seq[i], tag=mpi.anytag, comm=comm)
            tag <- mpi.get.sourcetag()[2]
            out[[tag]]<- tmp
            j <- i+slave.num
            if (j <= n)
                mpi.send.Robj(list(data.arg=list(x[[j]])), dest=apply.seq[i],tag=j, comm=comm)
            else
                mpi.send.Robj(as.integer(0),dest=apply.seq[i],tag=j,comm=comm)
        }
        return(out)
    }
    #.mpi.applyLB <<- integer(n)
	mpi.seq.tmp <- integer(n)
    for (i in 1:n){
       tmp<- mpi.recv.Robj(mpi.anysource,mpi.anytag,comm)
       srctag <- mpi.get.sourcetag()
       out[[srctag[2]]]<- tmp
       mpi.seq.tmp[i] <- srctag[1]
       j <- i+slave.num
       if (j <= n)
            mpi.send.Robj(list(data.arg=list(x[[j]])), dest=srctag[1],tag=j, comm=comm)
       else
            mpi.send.Robj(as.integer(0),dest=srctag[1],tag=j,comm=comm)
    }
 	#assign(".mpi.applyLB",mpi.seq.tmp, envir = .GlobalEnv)
	out
}

.mpi.worker.applyLB <- function(n){
    #assign(".mpi.err", FALSE,  envir = .GlobalEnv)
	.comm <- 1
    #n <- mpi.bcast(integer(1), type=1, comm=.comm)
    tmpfunarg <- mpi.bcast.Robj(rank=0, comm=.comm)
    .tmpfun <- tmpfunarg$fun
    dotarg <- tmpfunarg$dot.arg
    mpi.anytag <- mpi.any.tag()
    repeat {
        tmpdata.arg <- mpi.recv.Robj(source=0,tag=mpi.anytag, comm=.comm)$data.arg
        tag <- mpi.get.sourcetag()[2]
        if (tag > n)
            break
        out <- try(do.call(".tmpfun", c(tmpdata.arg, dotarg)),TRUE)
        #if (.mpi.err)
        #    print(geterrmessage())
        mpi.send.Robj(out,0,tag,.comm)
    }
}

mpi.iapplyLB <- function(x, fun, ...,  apply.seq=NULL, comm=1, sleep=0.01){
	apply.seq=NULL
    n <- length(x)
    slave.num <- mpi.comm.size(comm)-1
    if (slave.num < 1)
        stop("There are no slaves running")
    if (n <= slave.num) {
        if (exists(".mpi.applyLB"))
            rm(".mpi.applyLB",  envir =.GlobalEnv)
        return (mpi.iapply(x,fun,...,comm=comm,sleep=sleep))
    }
    if (!is.function(fun))
        stop("fun is not a function")
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
    mpi.bcast.Robj(list(fun=fun,dot.arg=list(...)),rank=0,comm=comm)
    out <- as.list(integer(n))
    mpi.anysource <- mpi.any.source()
    mpi.anytag <- mpi.any.tag()
    for (i in 1:slave.num)
        mpi.send.Robj(list(data.arg=list(x[[i]])), dest=i,tag=i,comm=comm)
    #for (i in 1:slave.num)
    #    mpi.waitany(slave.num)

    if (!is.null(apply.seq)){
       i=0
       repeat {
        if (mpi.iprobe(apply.seq[i+1],mpi.anytag,comm)){
            i=i+1
            j <- i+slave.num
            if ( j <= n)
                mpi.send.Robj(list(data.arg=list(x[[j]])), dest=apply.seq[i],tag=j, comm=comm) 
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
                mpi.send.Robj(list(data.arg=list(x[[j]])), dest=src,tag=j, comm=comm)
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
#        out <- try(do.call(".tmpfun", c(tmpdata.arg, dotarg)),TRUE)
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

mpi.parLapply <- function(x, fun, ..., job.num=mpi.comm.size(comm)-1, apply.seq=NULL, comm=1){
    if (job.num < 2)
        stop("job.num is at least 2.")
    splitList <- function(x, ncl)
        lapply(.splitIndices(length(x), ncl), function(i) x[i])
    .docall(c, mpi.applyLB(splitList(x, job.num), 
        lapply, fun, ..., apply.seq=apply.seq, comm=comm))
}

mpi.iparLapply <- function(x, fun, ..., job.num=mpi.comm.size(comm)-1, apply.seq=NULL, comm=1, sleep=0.01){
    if (job.num < 2)
        stop("job.num is at least 2.")
    splitList <- function(x, ncl)
        lapply(.splitIndices(length(x), ncl), function(i) x[i])
    .docall(c, mpi.iapplyLB(splitList(x, job.num), 
        lapply, fun, ..., apply.seq=apply.seq, comm=comm, sleep=sleep))
}

mpi.parSapply <- function (x, fun, ..., job.num=mpi.comm.size(comm)-1, apply.seq=NULL, 
                simplify = TRUE, USE.NAMES = TRUE, comm=1) 
{
    FUN <- match.fun(fun)
    answer <- mpi.parLapply(as.list(x),FUN,...,job.num=job.num,apply.seq=apply.seq,comm=comm)
    if (USE.NAMES && is.character(x) && is.null(names(answer))) 
        names(answer) <- x
    .simplify(length(x),answer, simplify)
}

mpi.iparSapply <- function (x, fun, ..., job.num=mpi.comm.size(comm)-1, apply.seq=NULL, 
                simplify = TRUE, USE.NAMES = TRUE, comm=1,sleep=0.01) 
{
    FUN <- match.fun(fun)
    answer <- mpi.iparLapply(as.list(x),FUN,...,job.num=job.num,apply.seq=apply.seq,comm=comm,sleep=sleep)
    if (USE.NAMES && is.character(x) && is.null(names(answer))) 
        names(answer) <- x
    .simplify(length(x),answer, simplify)
}

mpi.parReplicate <- function(n,  expr, job.num=mpi.comm.size(comm)-1, apply.seq=NULL,
                                simplify = TRUE, comm=1){
    mpi.parSapply(integer(n), eval.parent(substitute(function(...) expr)), 
    job.num=job.num, apply.seq=apply.seq, simplify = simplify, comm=comm)
}

mpi.iparReplicate <- function(n,  expr, job.num=mpi.comm.size(comm)-1, apply.seq=NULL,
                                simplify = TRUE, comm=1,sleep=0.01){
    mpi.iparSapply(integer(n), eval.parent(substitute(function(...) expr)), 
    job.num=job.num, apply.seq=apply.seq, simplify = simplify, comm=comm,sleep=sleep)
}

mpi.parRapply <- function(x,fun,...,job.num=mpi.comm.size(comm)-1,apply.seq=NULL,comm=1){
    if (job.num < 2)
        stop("job.num is at least 2.")
    splitRows <- function(x, ncl)
        lapply(.splitIndices(nrow(x), ncl), function(i) x[i,, drop=FALSE])
    .docall(c, mpi.applyLB(splitRows(x,job.num), apply, 1, fun, ..., 
            apply.seq=apply.seq, comm=comm))
}

mpi.iparRapply <- function(x,fun,...,job.num=mpi.comm.size(comm)-1,apply.seq=NULL,comm=1,sleep=0.01){
    if (job.num < 2)
        stop("job.num is at least 2.")
    splitRows <- function(x, ncl)
        lapply(.splitIndices(nrow(x), ncl), function(i) x[i,, drop=FALSE])
    .docall(c, mpi.iapplyLB(splitRows(x,job.num), apply, 1, fun, ..., 
            apply.seq=apply.seq, comm=comm,sleep=sleep))
}

mpi.parCapply <- function(x,fun,...,job.num=mpi.comm.size(comm)-1,apply.seq=NULL,comm=1){
    if (job.num < 2)
        stop("job.num is at least 2.")
    splitCols <- function(x, ncl)
        lapply(.splitIndices(ncol(x), ncl), function(i) x[,i, drop=FALSE])
    .docall(c, mpi.applyLB(splitCols(x,job.num), apply, 2, fun, ..., 
            apply.seq=apply.seq, comm=comm))
}

mpi.iparCapply <- function(x,fun,...,job.num=mpi.comm.size(comm)-1,apply.seq=NULL,comm=1,sleep=0.01){
    if (job.num < 2)
        stop("job.num is at least 2.")
    splitCols <- function(x, ncl)
        lapply(.splitIndices(ncol(x), ncl), function(i) x[,i, drop=FALSE])
    .docall(c, mpi.iapplyLB(splitCols(x,job.num), apply, 2, fun, ..., 
            apply.seq=apply.seq, comm=comm,sleep=sleep))
}

mpi.parApply <- function(x, MARGIN, fun, ..., job.num = mpi.comm.size(comm)-1,
                    apply.seq=NULL, comm=1)
{
    FUN <- match.fun(fun)
    d <- dim(x)
    dl <- length(d)
    if(dl == 0)
    stop("dim(x) must have a positive length")
    ds <- 1:dl

    if(length(oldClass(x)) > 0)
    x <- if(dl == 2) as.matrix(x) else as.array(x)
    dn <- dimnames(x)

    s.call <- ds[-MARGIN]
    s.ans  <- ds[MARGIN]
    d.call <- d[-MARGIN]
    d.ans  <- d[MARGIN]
    dn.call<- dn[-MARGIN]
    dn.ans <- dn[MARGIN]

    d2 <- prod(d.ans)
    if(d2 == 0) {
        newX <- array(vector(typeof(x), 1), dim = c(prod(d.call), 1))
        ans <- FUN(if(length(d.call) < 2) newX[,1] else
                   array(newX[,1], d.call, dn.call), ...)
        return(if(is.null(ans)) ans else if(length(d.call) < 2) ans[1][-1]
               else array(ans, d.ans, dn.ans))
    }
    newX <- aperm(x, c(s.call, s.ans))
    dim(newX) <- c(prod(d.call), d2)
    if(length(d.call) < 2) {
        if (length(dn.call)) dimnames(newX) <- c(dn.call, list(NULL))
        ans <- mpi.parLapply(1:d2, function(i, ...) FUN(newX[,i], ...), 
                ..., job.num = job.num, apply.seq=apply.seq, comm=comm )
    } 
    else
        ans <- mpi.parLapply(1:d2,
            function(i, ...) FUN(array(newX[,i], d.call, dn.call), ...),
                ..., job.num = job.num, apply.seq=apply.seq, comm=comm)
                    
    ans.list <- is.recursive(ans[[1]])
    l.ans <- length(ans[[1]])

    ans.names <- names(ans[[1]])
    if(!ans.list)
        ans.list <- any(unlist(lapply(ans, length)) != l.ans)
    if(!ans.list && length(ans.names)) {
        all.same <- sapply(ans, function(x) identical(names(x), ans.names))
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

mpi.iparApply <- function(x, MARGIN, fun, ..., job.num = mpi.comm.size(comm)-1,
                    apply.seq=NULL, comm=1,sleep=0.01)
{
    FUN <- match.fun(fun)
    d <- dim(x)
    dl <- length(d)
    if(dl == 0)
    stop("dim(x) must have a positive length")
    ds <- 1:dl

    if(length(oldClass(x)) > 0)
    x <- if(dl == 2) as.matrix(x) else as.array(x)
    dn <- dimnames(x)

    s.call <- ds[-MARGIN]
    s.ans  <- ds[MARGIN]
    d.call <- d[-MARGIN]
    d.ans  <- d[MARGIN]
    dn.call<- dn[-MARGIN]
    dn.ans <- dn[MARGIN]

    d2 <- prod(d.ans)
    if(d2 == 0) {
        newX <- array(vector(typeof(x), 1), dim = c(prod(d.call), 1))
        ans <- FUN(if(length(d.call) < 2) newX[,1] else
                   array(newX[,1], d.call, dn.call), ...)
        return(if(is.null(ans)) ans else if(length(d.call) < 2) ans[1][-1]
               else array(ans, d.ans, dn.ans))
    }
    newX <- aperm(x, c(s.call, s.ans))
    dim(newX) <- c(prod(d.call), d2)
    if(length(d.call) < 2) {
        if (length(dn.call)) dimnames(newX) <- c(dn.call, list(NULL))
        ans <- mpi.iparLapply(1:d2, function(i, ...) FUN(newX[,i], ...), 
                ..., job.num = job.num, apply.seq=apply.seq, comm=comm, sleep=sleep )
    } 
    else
        ans <- mpi.iparLapply(1:d2,
            function(i, ...) FUN(array(newX[,i], d.call, dn.call), ...),
                ..., job.num = job.num, apply.seq=apply.seq, comm=comm, sleep=sleep)
                    
    ans.list <- is.recursive(ans[[1]])
    l.ans <- length(ans[[1]])

    ans.names <- names(ans[[1]])
    if(!ans.list)
        ans.list <- any(unlist(lapply(ans, length)) != l.ans)
    if(!ans.list && length(ans.names)) {
        all.same <- sapply(ans, function(x) identical(names(x), ans.names))
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
