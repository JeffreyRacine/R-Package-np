test_that("density and distribution hats validate numeric RHS before coercion", {
  if (!spawn_mpi_slaves(1L)) skip("MPI slaves unavailable"); on.exit(close_mpi_slaves(),add=TRUE)
  withr::local_options(np.messages=FALSE)
  x <- data.frame(x=seq(.05,.95,length.out=24))
  ty <- data.frame(y=sin(seq_len(24)))
  fixtures <- list(
    npudenshat=list(bws=npudensbw(dat=x,bws=.3,bandwidth.compute=FALSE),tdat=x),
    npudisthat=list(bws=npudistbw(dat=x,bws=.3,bandwidth.compute=FALSE),tdat=x),
    npcdenshat=list(bws=npcdensbw(xdat=x,ydat=ty,bws=c(.3,.4),bandwidth.compute=FALSE),txdat=x,tydat=ty),
    npcdisthat=list(bws=npcdistbw(xdat=x,ydat=ty,bws=c(.3,.4),bandwidth.compute=FALSE),txdat=x,tydat=ty))
  bad <- list(factor(rep(c("1","2"),12)), ordered(rep(c("1","2"),12)),
              as.character(1:24), rep(TRUE,24), as.complex(1:24), as.list(1:24),
              data.frame(a=1:24,b=rep("2",24)), matrix(as.character(1:48),24,2))
  good <- list(1:24, as.double(1:24), matrix(1:48,24,2),
               data.frame(a=1:24,b=as.double(24:1)))
  for (owner in names(fixtures)) {
    fn <- getFromNamespace(owner, 'npRmpi')
    args <- fixtures[[owner]]
    H <- do.call(fn,c(args,list(output="matrix")))
    for (rhs in bad) for (output in c("matrix","apply"))
      expect_error(do.call(fn,c(args,list(y=rhs,output=output))),
                   "must have numeric (integer or double) storage",fixed=TRUE,info=owner)
    for (rhs in good) {
      expected <- H %*% as.matrix(rhs)
      actual <- do.call(fn,c(args,list(y=rhs,output="apply")))
      expect_equal(as.numeric(actual),as.numeric(expected),tolerance=2e-11,info=owner)
      attached <- do.call(fn,c(args,list(y=rhs,output="matrix")))
      expect_equal(as.numeric(attr(attached,"Hy")),as.numeric(expected),tolerance=2e-11,info=owner)
    }
    rhs <- as.double(1:24); rhs[3L] <- NA_real_
    args.small <- args
    for (name in intersect(c("tdat","txdat","tydat"),names(args.small)))
      args.small[[name]] <- args.small[[name]][-3L,,drop=FALSE]
    small <- do.call(fn,c(args.small,list(y=rhs[-3L],output="apply")))
    with.na <- do.call(fn,c(args,list(y=rhs,output="apply")))
    expect_equal(as.numeric(with.na),as.numeric(small),tolerance=2e-11,info=owner)
  }
})
