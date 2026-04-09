.npRmpi_validate_raw_length <- function(n, caller) {
    if (length(n) != 1L || is.na(n))
        stop(sprintf("%s: expected a single non-missing length value", caller))
    n <- as.integer(n)
    if (n < 0L)
        stop(sprintf("%s: expected a non-negative length value", caller))
    n
}

.npRmpi_validate_rcounts <- function(rcounts, size, total_len, caller) {
    if (length(rcounts) != size)
        stop(sprintf("%s: expected %d receive counts, got %d", caller, size, length(rcounts)))
    if (any(is.na(rcounts)))
        stop(sprintf("%s: receive counts contain NA", caller))
    rcounts <- as.integer(rcounts)
    if (any(rcounts < 0L))
        stop(sprintf("%s: receive counts must be non-negative", caller))
    if (sum(rcounts) != total_len)
        stop(sprintf("%s: receive counts sum %d but payload length is %d",
                     caller, sum(rcounts), total_len))
    rcounts
}

.npRmpi_split_raw_by_counts <- function(allbiobj, rcounts, caller) {
    size <- length(rcounts)
    rcounts <- .npRmpi_validate_rcounts(rcounts, size, length(allbiobj), caller)
    pos <- c(0L, cumsum(rcounts))
    cutobj <- vector("list", size)
    for (i in seq_len(size)) {
        if (rcounts[i] == 0L) {
            cutobj[[i]] <- raw(0)
        } else {
            start <- pos[i] + 1L
            stopifnot(start <= pos[i + 1L])
            cutobj[[i]] <- allbiobj[start:pos[i + 1L]]
        }
    }
    cutobj
}

### Copyright (C) 2002 Hao Yu
mpi.probe <- function(source, tag, comm=1, status=0){
    .Call("mpi_probe", as.integer(source), as.integer(tag), 
        as.integer(comm), as.integer(status),
        PACKAGE = "npRmpi")
}

mpi.get.count <- function(type, status = 0){
    .Call("mpi_get_count",as.integer(status), 
        as.integer(type),PACKAGE = "npRmpi")
}

mpi.get.sourcetag <- function(status=0){
    .Call("mpi_get_sourcetag", as.integer(status),PACKAGE = "npRmpi")
}

mpi.gather <- function(x, type, rdata, root=0,  comm=1){
    .Call("mpi_gather",.force.type( x,type), as.integer(type), rdata, as.integer(root), as.integer(comm),PACKAGE = "npRmpi")
}

mpi.scatter <- function(x, type, rdata, root=0,  comm=1){
    .Call("mpi_scatter", .force.type(x,type), as.integer(type), rdata, 
        as.integer(root), as.integer(comm),PACKAGE = "npRmpi")
}

mpi.gatherv <- function(x, type, rdata, rcounts, root=0,  comm=1){
    .Call("mpi_gatherv", x, as.integer(type),rdata, as.integer(rcounts), 
        as.integer(root), as.integer(comm),PACKAGE = "npRmpi")
}

mpi.scatterv <- function(x, scounts, type, rdata, root=0, comm=1){
    .Call("mpi_scatterv", .force.type(x,type), as.integer(scounts), as.integer(type), rdata, 
        as.integer(root), as.integer(comm),PACKAGE = "npRmpi")
}
#strings.link<-function(manysts,newst){
#.Call("stringslink",as.character(manysts),as.character(newst),PACKAGE="npRmpi")
#}

#string.cut<-function(onest,newst){
#.Call("stringcut",as.character(onest),as.character(newst),PACKAGE="npRmpi")
#}



mpi.scatter.Robj <- function(obj=NULL, root=0, comm=1){
    if (mpi.comm.rank(comm) == root){
			size<-mpi.comm.size(comm)
        #subobj<-lapply(obj,serialize, connection=NULL)
			subobj<-lapply(seq_len(size), function(i) serialize(obj[[i]], NULL))
	
			sublen<-unlist(lapply(subobj,length))
        #newsubobj<-strings.link(subobj,string(sum(sublen)+1))
			newsubobj<-c(subobj,recursive=TRUE)
        strnum <- .npRmpi_validate_raw_length(
            mpi.scatter(sublen,type=1,rdata=integer(1),root=root,comm=comm),
            "mpi.scatter.Robj")
			outobj<-unserialize(mpi.scatterv(newsubobj,scounts=sublen,type= 4,
                rdata=raw(strnum),root=root, comm=comm))
    }
    else {
        strnum <- .npRmpi_validate_raw_length(
            mpi.scatter(integer(1),type=1,rdata=integer(1),root=root,comm=comm),
            "mpi.scatter.Robj")
        outobj<-unserialize(mpi.scatterv(raw(strnum),scounts=0, type=4,
                rdata=raw(strnum), root=root, comm=comm))
    }
	gc()
    return(outobj)
}

mpi.scatter.Robj2slave=function (obj, comm = 1) {
    if (!is.list(obj))
        stop("Only list object is allowed to scatter to slaves.")
    if (length(obj) != (mpi.comm.size(comm)-1)) 
        stop("Length of your list object is not the same as total number of slaves.")
    .tmpname <- list(objname=deparse(substitute(obj), width.cutoff = 500))
    mpi.bcast.Robj2slave(.tmpname)
    mpi.bcast.cmd(cmd = .tmpRobj <- mpi.scatter.Robj(comm = 1), 
        rank = 0, comm = comm)
    mpi.scatter.Robj(obj=c(list("master"),obj), root = 0, comm = comm)
    mpi.bcast.cmd(cmd = assign(.tmpname$objname, .tmpRobj), rank = 0, comm = comm)
}

mpi.gather.Robj <- function(obj=NULL, root=0, comm=1, ...){
    biobj<-serialize(obj, NULL)
    bilen<-length(biobj)
    if (mpi.comm.rank(comm) == root){
        size<-mpi.comm.size(comm=comm)
        rcounts<-mpi.gather(bilen,type=1,rdata=integer(size),
            root=root,comm=comm)
        allbiobj<-mpi.gatherv(biobj,type=4,rdata=raw(sum(rcounts))
                        ,rcounts=rcounts,root=root,comm=comm)
    cutobj <- .npRmpi_split_raw_by_counts(allbiobj, rcounts, "mpi.gather.Robj")
			out <- sapply(cutobj,unserialize, ..., USE.NAMES = FALSE)
			gc()
			out
    }
    else {
         mpi.gather(bilen,type=1,rdata=integer(1),root=root,comm=comm)
         out <- mpi.gatherv(biobj,type=4,rdata=raw(1),rcounts=0,root=root,comm=comm)
		 gc()
		 out
   }
}

mpi.allgather.Robj <- function(obj=NULL, comm=1){
    biobj<-serialize(obj, NULL)
    bilen<-length(biobj)
    size<-mpi.comm.size(comm=comm)
    rcounts<-mpi.allgather(bilen,type=1,rdata=integer(size),comm=comm)
    allbiobj<-mpi.allgatherv(biobj,type=4,rdata=raw(sum(rcounts))
        ,rcounts=rcounts,comm=comm)
    cutobj <- .npRmpi_split_raw_by_counts(allbiobj, rcounts, "mpi.allgather.Robj")
    out <- sapply(cutobj,unserialize, USE.NAMES = FALSE)
	gc()
	out
   # bistrcut<-sapply(rcounts,string)
   # bistr<-string.cut(allbiobj,bistrcut)
   # lapply(bistr,unserialize)
}



mpi.allgather <- function(x, type, rdata, comm=1){
    .Call("mpi_allgather", x, as.integer(type), rdata, as.integer(comm),
        PACKAGE = "npRmpi")
}

mpi.allgatherv <- function(x, type, rdata, rcounts, comm=1){
    .Call("mpi_allgatherv", x, as.integer(type), rdata, 
    as.integer(rcounts), as.integer(comm),PACKAGE = "npRmpi")
}

mpi.bcast <- function (x, type, rank = 0, comm = 1, buffunit=100) {
    .Call("mpi_bcast", .force.type(x,type), as.integer(type), as.integer(rank), 
        as.integer(comm), as.integer(buffunit), PACKAGE = "npRmpi")
}

.npRmpi_bcast_cmd_ns_funref <- function(scmd) {
    if (!is.call(scmd) || length(scmd) < 3L)
        return(NULL)
    hd <- scmd[[1L]]
    if (!is.symbol(hd) || !as.character(hd) %in% c("::", ":::"))
        return(NULL)

    pkg <- scmd[[2L]]
    tgt <- scmd[[3L]]
    pkg_name <- if (is.symbol(pkg)) as.character(pkg) else if (is.character(pkg) && length(pkg) >= 1L) pkg[[1L]] else ""
    tgt_name <- if (is.symbol(tgt)) as.character(tgt) else if (is.character(tgt) && length(tgt) >= 1L) tgt[[1L]] else ""
    if (!nzchar(pkg_name) || !nzchar(tgt_name))
        return(NULL)

    if (!requireNamespace(pkg_name, quietly = TRUE))
        return(NULL)

    if (as.character(hd) == "::")
        return(tryCatch(getExportedValue(pkg_name, tgt_name), error = function(e) NULL))

    tryCatch(get(tgt_name, envir = asNamespace(pkg_name), mode = "function", inherits = FALSE),
             error = function(e) NULL)
}

.npRmpi_bcast_cmd_lookup_fun <- function(name, eval_env = parent.frame()) {
    if (!is.character(name) || length(name) < 1L || !nzchar(name[[1L]]))
        return(NULL)
    nm <- name[[1L]]
    fn <- tryCatch(get(nm, mode = "function", envir = eval_env, inherits = TRUE),
                   error = function(e) NULL)
    if (is.function(fn))
        return(fn)

    tryCatch(get(nm, mode = "function", envir = asNamespace("npRmpi"), inherits = FALSE),
             error = function(e) NULL)
}

.npRmpi_bcast_cmd_funref <- function(scmd, eval_env = parent.frame()) {
    if (is.function(scmd))
        return(scmd)
    if (is.symbol(scmd)) {
        nm <- as.character(scmd)
        fn <- .npRmpi_bcast_cmd_lookup_fun(nm, eval_env = eval_env)
        if (is.function(fn))
            return(fn)
        return(nm)
    }
    if (is.character(scmd) && length(scmd) >= 1L) {
        fn <- .npRmpi_bcast_cmd_lookup_fun(scmd[[1L]], eval_env = eval_env)
        if (is.function(fn))
            return(fn)
        return(scmd[[1L]])
    }
    if (is.call(scmd) && length(scmd) >= 1L) {
        hd <- scmd[[1L]]
        if (is.symbol(hd) && as.character(hd) %in% c("::", ":::") && length(scmd) >= 3L) {
            fn <- .npRmpi_bcast_cmd_ns_funref(scmd)
            if (is.function(fn))
                return(fn)
        }
        if (is.call(hd)) {
            fn <- tryCatch(.npRmpi_eval_scmd(hd, envir = eval_env), error = function(e) NULL)
            if (is.function(fn))
                return(fn)
        }
        if (is.function(hd))
            return(hd)
        if (is.symbol(hd)) {
            nm <- as.character(hd)
            fn <- .npRmpi_bcast_cmd_lookup_fun(nm, eval_env = eval_env)
            if (is.function(fn))
                return(fn)
            return(nm)
        }
        if (is.character(hd) && length(hd) >= 1L) {
            fn <- .npRmpi_bcast_cmd_lookup_fun(hd[[1L]], eval_env = eval_env)
            if (is.function(fn))
                return(fn)
            return(hd[[1L]])
        }
    }
    as.character(scmd)[1L]
}

.npRmpi_bcast_cmd_is_plot_call <- function(scmd) {
    if (!is.call(scmd) || length(scmd) < 1L)
        return(FALSE)
    hd <- scmd[[1L]]
    if (is.symbol(hd))
        return(identical(as.character(hd), "plot"))
    if (is.call(hd) && length(hd) >= 3L && is.symbol(hd[[1L]]) &&
        as.character(hd[[1L]]) %in% c("::", ":::")) {
        tgt <- hd[[3L]]
        return(is.symbol(tgt) && identical(as.character(tgt), "plot"))
    }
    FALSE
}

#bin.nchar <- function(x){
#    if (!is.character(x))
#        stop("Must be a (binary) character")
#    .Call("bin_nchar", x[1],PACKAGE = "npRmpi")
#}

mpi.bcast.cmd <- function (cmd=NULL, ..., rank=0, comm=1, nonblock=FALSE, sleep=0.1, caller.execute = FALSE){
	myrank=mpi.comm.rank(comm)
    if(myrank == rank){
      scmd <- substitute(cmd)
      if (isTRUE(.npRmpi_bcast_cmd_is_plot_call(scmd))) {
          stop("plot(...) inside mpi.bcast.cmd(...) is unsupported in canonical SPMD mode; call plot(...) directly from master with npRmpi.autodispatch=TRUE", call. = FALSE)
      }
      if(caller.execute) tcmd <- substitute(cmd)
        #cmd <- deparse(substitute(cmd), width.cutoff=500)
		#cmd <- serialize(cmd, NULL)
		arg <-list(...)
		commsize <- mpi.comm.size(comm=comm)

		scmd.arg <-serialize(list(scmd=scmd, arg=arg), NULL)
		#mpi.bcast(x=length(cmd), type=1, rank=rank, comm=comm)
		#invisible(mpi.bcast(x=cmd, type=4, rank=rank, comm=comm))
		
		for (i in 0:(commsize - 1)) {
			if (i != rank)
				invisible(mpi.send(x=scmd.arg, type=4, dest=i, tag=50000+i, comm=comm))
			}
	      if (caller.execute) {
          .npRmpi_with_manual_bcast_context({
	       if (length(arg) > 0)
	         do.call(.npRmpi_bcast_cmd_funref(scmd), arg, envir = parent.frame())
	       else
	         .npRmpi_eval_scmd(tcmd, envir = parent.frame())
          })
	      }
  
    } 
    else {
       # charlen <- mpi.bcast(x=integer(1), type=1, rank=rank, comm=comm)
		if (!nonblock){
			mpi.probe(mpi.any.source(), tag=50000+myrank, comm)
			srctag <- mpi.get.sourcetag(0)
			charlen <- mpi.get.count(type=4, 0)
			#out <- unserialize(mpi.recv(x=raw(charlen), type=4,srctag[1],srctag[2], comm))
			scmd.arg <- unserialize(mpi.recv(x=raw(charlen), type=4,srctag[1],srctag[2], comm))

		} else {
			repeat {
				if (mpi.iprobe(mpi.any.source(),tag=50000+myrank,comm)){
					srctag <- mpi.get.sourcetag()
					charlen <- mpi.get.count(type=4)
					#out <- unserialize(mpi.recv(x = raw(charlen), type = 4, srctag[1],srctag[2], comm))
					scmd.arg <- unserialize(mpi.recv(x = raw(charlen), type = 4, srctag[1],srctag[2], comm))
					break
				}
				Sys.sleep(sleep)
			}
		}
				if (length(scmd.arg$arg)>0) {
				as.call(list(
            .npRmpi_bcast_cmd_funref(".npRmpi_with_manual_bcast_context"),
				  as.call(list(
				    as.name("do.call"),
				    .npRmpi_bcast_cmd_funref(scmd.arg$scmd),
				    scmd.arg$arg,
				    quote = FALSE,
				    envir = .GlobalEnv
				  ))
				))
			} else
        as.call(list(
          .npRmpi_bcast_cmd_funref(".npRmpi_with_manual_bcast_context"),
          scmd.arg$scmd
        ))
		}
}

mpi.bcast.Robj <- function(obj=NULL, rank=0, comm=1){
    if (mpi.comm.rank(comm) == rank){
    tmp <- serialize(obj, NULL)
    mpi.bcast(as.integer(length(tmp)), 1, rank, comm)
    mpi.bcast(tmp, 4, rank, comm)
	invisible(NULL)
    }
    else {
    charlen <- mpi.bcast(integer(1), 1, rank, comm)
    out <- unserialize(mpi.bcast(raw(charlen), 4, rank, comm))
	out
    }
}

mpi.bcast.Robj2slave <- function(obj, comm=1, all=FALSE){
    if (!all){
		objname <- deparse(substitute(obj),width.cutoff=500)
		obj <- list(objname=objname,obj=obj)
		mpi.bcast.cmd(cmd=.tmpRobj <- mpi.bcast.Robj(comm=1),
                    rank=0, comm=comm)
		mpi.bcast.Robj(obj, rank=0, comm=comm)
		mpi.bcast.cmd(cmd=assign(.tmpRobj$objname,.tmpRobj$obj), rank=0, comm=comm)
		#mpi.bcast.cmd(rm(.tmpRobj,envir = .GlobalEnv), rank=0, comm=comm) 
	}
	else {
		master.objects <- objects(envir = .GlobalEnv)
		object.values <- mget(master.objects,
		                      envir = .GlobalEnv,
		                      inherits = FALSE,
		                      ifnotfound = vector("list", length(master.objects)))
		obj.num=length(master.objects)
		if (obj.num)
			for (i in seq_len(obj.num)){
				mpi.bcast.cmd(cmd=.tmpRobj <- mpi.bcast.Robj(comm=1),
	                    rank=0, comm=comm)
				mpi.bcast.Robj(list(objname=master.objects[i], obj=object.values[[i]]), 
					rank=0, comm=comm)
				mpi.bcast.cmd(cmd=assign(.tmpRobj$objname,.tmpRobj$obj), rank=0, comm=comm)
			}
	}
}

mpi.bcast.Rfun2slave <- function(comm=1){
	master.fun <- objects(envir = .GlobalEnv)
	fun.values <- mget(master.fun,
	                   envir = .GlobalEnv,
	                   inherits = FALSE,
	                   ifnotfound = vector("list", length(master.fun)))
	sync.index <- which(vapply(fun.values, is.function, logical(1)))
	obj.num=length(sync.index)
	if (obj.num)
		for (i in sync.index){
			mpi.bcast.cmd(cmd=.tmpRobj <- mpi.bcast.Robj(comm=1),
	                   rank=0, comm=comm)
			mpi.bcast.Robj(list(objname=master.fun[i], obj=fun.values[[i]]), 
				rank=0, comm=comm)
			mpi.bcast.cmd(cmd=assign(.tmpRobj$objname,.tmpRobj$obj), rank=0, comm=comm)
		}
}

mpi.bcast.data2slave <- function(obj, comm=1, buffunit=100){
	if (!is.numeric(obj) || (!is.vector(obj) && !is.matrix(obj)))
		return (mpi.bcast.Robj2slave(obj, comm=comm))
		#stop ("Please use mpi.bcast.Robj2slave")
		
	objname <- serialize(deparse(substitute(obj),width.cutoff=500),NULL)
	obj.info = integer(4)
	obj.info[1]=length(objname)
	if (is.vector(obj)){
		if (buffunit < 1 || buffunit >=2^31)
			stop("buffunit must be an integer between 1 and 2^31-1")
		obj.info[2]=buffunit
		obj.info[3]=length(obj)%/%buffunit
		obj.info[4]=length(obj)%%buffunit
	}
	if (is.matrix(obj)){
		obj.dim <-dim(obj)
		obj.info[2]=obj.dim[1]
		obj.info[3]=obj.dim[2]
		obj.info[4]=0
	}
	
	mpi.bcast.cmd(.tinfo <- mpi.bcast(integer(4),type=1),rank=0,comm=1)
	mpi.bcast(obj.info,type=1,rank=0,comm=comm)
	
	mpi.bcast.cmd(.tname<-unserialize(mpi.bcast(raw(.tinfo[1]),type=4)), rank=0, comm=comm)
	mpi.bcast(objname, type=4, rank=0, comm=comm)
	
	if (is.vector(obj)){
		mpi.bcast.cmd(.tmp.obj <- mpi.bcast(double(.tinfo[2]*(.tinfo[3]+(.tinfo[4]>0))),type=5, buffunit=.tinfo[2]),rank=0,comm=comm)
		mpi.bcast(obj,type=5,rank=0,comm=comm,buffunit=buffunit)
			mpi.bcast.cmd(assign(.tname,.tmp.obj[seq_len(.tinfo[2]*.tinfo[3]+.tinfo[4])]), rank=0, comm=comm)
		mpi.bcast.cmd(rm(".tmp.obj"))
	}
	if (is.matrix(obj)){
		mpi.bcast.cmd(.tmp.obj <- mpi.bcast(matrix(double(.tinfo[2]*.tinfo[3]), nrow=.tinfo[2]),type=5, buffunit=.tinfo[2]),rank=0,comm=comm)
		mpi.bcast(obj,type=5,rank=0,comm=comm,buffunit=obj.info[2])	
		mpi.bcast.cmd(assign(.tname,.tmp.obj), rank=0, comm=comm)
	}
}

mpi.send <- function (x, type,  dest, tag, comm=1){
    .Call("mpi_send", .force.type(x,type), as.integer(type), as.integer(dest), 
    as.integer(tag), as.integer(comm),PACKAGE = "npRmpi")
}

mpi.recv <- function (x, type, source, tag, comm=1, status=0){
    .Call("mpi_recv", x, as.integer(type), as.integer(source), 
    as.integer(tag), as.integer(comm), as.integer(status),
    PACKAGE = "npRmpi")
}

mpi.send.Robj <- function(obj, dest, tag, comm=1){
    mpi.send(x=serialize(obj, NULL), type=4, dest=dest, tag=tag, comm=comm)
	invisible(NULL)
}

mpi.recv.Robj <- function(source, tag, comm=1, status=0){
    mpi.probe(source, tag, comm, status)
    srctag <- mpi.get.sourcetag(status)
    charlen <- mpi.get.count(type=4, status)
    out<-unserialize(mpi.recv(x=raw(charlen), type=4,srctag[1],srctag[2], comm, status))
	#gc()
	out
}

mpi.reduce <- function(x, type=2, 
    op=c("sum","prod","max","min","maxloc","minloc"), dest=0, comm=1){
#   op <- switch(match.arg(op),sum=1,prod=2,max=3,min=4)
    op <- pmatch(match.arg(op), 
        c("sum","prod","max","min","maxloc","minloc"))
    if (is.integer(x)){
       if(type!=1)
        stop("data (integer) and type are not matched.")
    }
    else if (is.double(x)){
       if(type!=2)
        stop("data (double) and type are not matched.")
    }
    else 
        stop("Not implemented.")

#      if (op==5||op==6){
#           n <- length(x)
#           x <- rep(x,rep(2,n))
#       x[seq(2, 2*n, 2)] <- mpi.comm.rank(comm)
#   }
        
    .Call("mpi_reduce", x, as.integer(type), as.integer(op), 
        as.integer(dest), as.integer(comm),PACKAGE = "npRmpi")
}

mpi.allreduce <- function(x,type=2,
    op=c("sum","prod","max","min","maxloc","minloc"), comm=1){
#   op <- switch(match.arg(op),sum=1,prod=2,max=3,min=4)
    op <- pmatch(match.arg(op), c("sum","prod","max","min","maxloc","minloc"))
    if (is.integer(x)){
       if(type!=1)
        stop("data (integer) and type are not matched.")
    }
    else if (is.double(x)){
       if(type!=2)
        stop("data (double) and type are not matched.")
    }
    else 
        stop("Not implemented.")
    .Call("mpi_allreduce", x, as.integer(type), as.integer(op), 
        as.integer(comm),PACKAGE = "npRmpi")
}

mpi.isend <- function (x, type,  dest, tag, comm=1, request=0){
    stop("mpi.isend is temporarily unsupported in npRmpi; use blocking mpi.send() or mpi.send.Robj() instead")
}

mpi.irecv <- function (x, type, source, tag, comm=1, request=0){
    stop("mpi.irecv is temporarily unsupported in npRmpi; use blocking mpi.recv() or mpi.recv.Robj() instead")
}

mpi.isend.Robj <- function(obj, dest, tag, comm=1,request=0){
    stop("mpi.isend.Robj is temporarily unsupported in npRmpi; use blocking mpi.send.Robj() instead")
}

mpi.wait <- function(request, status=0)
    invisible(.Call("mpi_wait",  as.integer(request), as.integer(status), PACKAGE = "npRmpi"))
    
mpi.waitany <- function(count, status=0){
    #mpi.realloc.request(count)
    .Call("mpi_waitany",  as.integer(count), as.integer(status), PACKAGE = "npRmpi")
}

mpi.waitall <- function(count){
    #mpi.realloc.request(count)
    #mpi.realloc.status(count)
    invisible(.Call("mpi_waitall",  as.integer(count), PACKAGE = "npRmpi"))
}

mpi.waitsome <- function(count){
    #mpi.realloc.request(count)
    #mpi.realloc.status(count)
    tmp<-.Call("mpi_waitsome",  as.integer(count), PACKAGE = "npRmpi")
    if (tmp[1] <0 || tmp[1] > count)
        return(list(count=tmp[1],indices=NULL))
    else 
        return(list(count=tmp[1],indices=tmp[2:(1+tmp[1])]))
}

mpi.test <- function(request, status=0)
    as.logical(.Call("mpi_test",  as.integer(request), as.integer(status), PACKAGE = "npRmpi"))

mpi.testany <- function(count, status=0){
    #mpi.realloc.request(count)
    tmp <-.Call("mpi_testany",  as.integer(count), as.integer(status), PACKAGE = "npRmpi")
    list(index=tmp[1], flag=as.logical(tmp[2]))
}

mpi.testall <- function(count){
    #mpi.realloc.request(count)
    #mpi.realloc.status(count)
    as.logical(.Call("mpi_testall",  as.integer(count), PACKAGE = "npRmpi"))
}

mpi.testsome <- function(count){
    #mpi.realloc.request(count)
    #mpi.realloc.status(count)
    tmp<-.Call("mpi_testsome",  as.integer(count), PACKAGE = "npRmpi")
    if (tmp[1] < 0 || tmp[1] > count)
        return(list(count=tmp[1],indices=NULL))
    else 
        return(list(count=tmp[1],indices=tmp[2:(1+tmp[1])]))
}

mpi.cancel <- function(request)
    invisible(.Call("mpi_cancel",  as.integer(request), PACKAGE = "npRmpi"))

mpi.test.cancelled <- function(status=0)
    as.logical(.Call("mpi_test_cancelled", as.integer(status), PACKAGE = "npRmpi"))

mpi.iprobe <- function(source, tag, comm=1, status=0){
    as.logical(.Call("mpi_iprobe", as.integer(source), as.integer(tag), 
        as.integer(comm), as.integer(status),
        PACKAGE = "npRmpi"))
}

mpi.realloc.status <- function(newmaxsize)
    if (newmaxsize > mpi.status.maxsize())
        invisible(.Call("mpi_realloc_status", as.integer(newmaxsize), PACKAGE = "npRmpi"))

mpi.realloc.request <- function(newmaxsize)
    if (newmaxsize > mpi.request.maxsize())
        invisible(.Call("mpi_realloc_request", as.integer(newmaxsize), PACKAGE = "npRmpi"))

mpi.realloc.comm <- function(newmaxsize)
    if (newmaxsize > mpi.comm.maxsize())
       invisible(.Call("mpi_realloc_comm", as.integer(newmaxsize), PACKAGE = "npRmpi"))

mpi.comm.maxsize <- function()
    .Call("mpi_comm_maxsize", PACKAGE = "npRmpi")

mpi.status.maxsize <- function()
    .Call("mpi_status_maxsize", PACKAGE = "npRmpi")
    
mpi.request.maxsize <- function()
    .Call("mpi_request_maxsize", PACKAGE = "npRmpi")
    
.mpi.undefined <- function()
    .Call("mpi_undefined", PACKAGE = "npRmpi")

.force.type <- function(x, type){ 
    if (length(type) != 1L || is.na(type))
        stop(".force.type: 'type' must be a single non-missing integer code")
    type <- as.integer(type)
    if (!(type %in% 1:5))
        return(x)
    switch(type,
        as.integer(x),
        as.double(x),
        as.character(x),
			as.raw(x),
		as.double(x))
}
#.mpi.serialize<- function(obj){
#    trans_obj=serialize(obj,NULL)
#    if ( version$major > 2 || version$minor >= 4.0)
#    if (getRversion()>="2.4.0")
#        return(trans_obj)
#    else
#        return(charToRaw(trans_obj))
#}


#.mpi.unserialize<- function(obj){
#    #if ( version$major > 2 || version$minor >= 4.0)
#    if (getRversion()>="2.4.0")
#        return(unserialize(obj))
#    else
#        return(unserialize(rawToChar(obj)))
# }
#mpi.request.get.status <- function(request, status=0){
#    as.logical(.Call("mpi_request_get_status",  as.integer(request), 
#        as.integer(status), PACKAGE = "npRmpi"))
#}
