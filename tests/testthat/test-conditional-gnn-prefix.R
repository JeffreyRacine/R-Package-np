cgnn_fixture <- function() list(
  x=c(-.67932861763984,.597772665787488,-.135870280675590,-.180028952658176,
      -.0812444845214486,-.0922158942557871,-.119163474533707,.565201563760638,
      -.167390177026391,-.617473398335278,-.452014316339046,.213771468494087),
  y=c(.117393371024781,-.215783902381389,-.752371990755028,1.55294725882662,
      .776713307792541,.329678953352629,-1.89405512568533,-1.99349412668255,
      -.220072074888336,.198501257716297,.563714301921436,-.880794421005045))

# Independent literal deleted-sample WLS and full-line numerical integration.
# No package NN decoder, cached basis or influence helper is used here.
cgnn_literal <- function(x,y,kx,ky,degree,xkernel='epanechnikov',xorder=2L,yorder=2L,ykernel='gaussian') {
  response <- function(z) {
    v <- z*z
    if(ykernel=='uniform')return(.5*(v<1))
    if(ykernel=='epanechnikov'){
      out<-switch(as.character(yorder),
        '2'=.33541019662496845446*(1-v/5),
        '4'=.008385254916*(-15+7*v)*(-5+v),
        '6'=.33541019662496845446*(2.734375+v*(-3.28125+.721875*v))*(1-.2*v),
        '8'=.33541019662496845446*(3.5888671875+v*(-7.8955078125+
          v*(4.1056640625-.5865234375*v)))*(1-.2*v))
      out[v>=5]<-0;return(out)
    }
    polynomial <- switch(as.character(yorder), '2'=rep(1,length(z)),
      '4'=1.5-.5*v, '6'=1.875+v*(-1.25+.125*v),
      '8'=2.1875+v*(-2.1875+v*(.4375-.02083333333*v)))
    stats::dnorm(z)*polynomial
  }
  n <- length(y); terms <- 0:degree
  basis <- outer(x,terms,`^`)
  rows <- lapply(seq_len(n),function(i) {
    keep <- seq_len(n)!=i
    lookup <- min(kx,n-2L)
    h <- sort(abs(x[i]-x[keep]))[lookup]*(kx/lookup)
    w <- if(xkernel=='uniform') .5*(abs(x[i]-x[keep])<h) else
      3/(4*sqrt(5))*pmax(0,1-((x[i]-x[keep])/h)^2/5)
    if(xkernel=='epanechnikov' && xorder!=2L) {
      v<-((x[i]-x[keep])/h)^2
      w<-switch(as.character(xorder),
        '4'=.008385254916*(-15+7*v)*(-5+v),
        '6'=.33541019662496845446*(2.734375+v*(-3.28125+.721875*v))*(1-.2*v),
        '8'=.33541019662496845446*(3.5888671875+v*(-7.8955078125+
          v*(4.1056640625-.5865234375*v)))*(1-.2*v))
      w[v>=5]<-0
    }
    B <- basis[keep,,drop=FALSE]
    w*as.vector(B%*%solve(crossprod(B,w*B),basis[i,]))
  })
  I2 <- mean(vapply(seq_len(n),function(i) {
    lookup <- min(ky,n-2L)
    h <- sort(abs(y[i]-y[-i]))[lookup]*(ky/lookup)
    sum(rows[[i]]*response((y[i]-y[-i])/h)/h)
  },0))
  cuts <- c(-Inf,sort(unique(c(y,as.vector(outer(y,y,`+`)/2)))),Inf)
  if(ykernel!='gaussian'){
    support<-(if(ykernel=='uniform')1 else sqrt(5))*ky/min(ky,n-2L)
    # Independent q-space anchor/donor equalities, not native reciprocal cuts.
    for(sign in c(-1,1))if(1-sign*support!=0)
      cuts<-c(cuts,as.vector(outer(y,sign*support*y,`-`)/(1-sign*support)))
    cuts<-sort(unique(cuts))
  }
  I1 <- sum(vapply(seq_len(length(cuts)-1L),function(j)
    integrate(function(q)vapply(q,function(t)
      mean(vapply(seq_len(n),function(i) {
        lookup <- min(ky,n-2L)
        h <- sort(abs(t-y[-i]))[lookup]*(ky/lookup)
        sum(rows[[i]]*response((t-y[-i])/h)/h)^2
      },0)),0),cuts[j],cuts[j+1L],abs.tol=1e-12/length(cuts),
      rel.tol=1e-10,subdivisions=1000L)$value,0))
  c(I1=I1,I2=I2,score=2*I2-I1)
}

test_that('conditional GNN local moments preserve small within-cluster distances', {
  old <- options(np.messages=FALSE,np.tree=FALSE);on.exit(options(old),add=TRUE)
  skip_if_not(spawn_mpi_slaves(1L), 'MPI session unavailable')
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  # Unique cancellation sentinel: ordinary globally scaled fixtures cannot
  # detect the loss of local X differences before compensated summation.
  # Larger n, all kernel orders and timing belong to the forensic campaign.
  x <- c(0,1e-8,2e-8,3e-8,1);y <- c(2,-1,1,0,-2)
  expected <- cgnn_literal(x,y,2L,3L,0L,ykernel='uniform')['score']
  for(tree in c(FALSE,TRUE)) {
    options(np.tree=tree)
    b <- npcdensbw(xdat=data.frame(x=x),ydat=y,bws=c(3,2),
      bwtype='generalized_nn',bwmethod='cv.ls',regtype='lc',
      cxkertype='epanechnikov',cykertype='uniform',bandwidth.compute=FALSE)
    command <- substitute({
      options(np.tree=TREE,np.messages=FALSE)
      get('.npcdensbw_eval_only',asNamespace('npRmpi'))(
        data.frame(x=X),Y,B,invalid.penalty='dbmax',force.local=FALSE)
    },list(TREE=tree,X=x,Y=y,B=b))
    observed <- get('.npRmpi_bcast_cmd_expr',asNamespace('npRmpi'))(
      command,comm=1L,caller.execute=TRUE)$objective
    expect_lte(abs(observed-expected),1e-9)
  }
})

test_that('conditional GNN prefix uses the literal whole-support criterion', {
  old <- options(np.messages=FALSE,np.largeh=TRUE,np.largelambda=TRUE)
  on.exit(options(old),add=TRUE)
  skip_if_not(spawn_mpi_slaves(1L), "MPI session unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  d <- cgnn_fixture()
  ev <- get('.npcdensbw_eval_only',asNamespace('npRmpi'))
  ref <- lapply(0:3,function(degree)cgnn_literal(d$x,d$y,8L,3L,degree))
  for(degree in 0:3) for(bernstein in c(FALSE,TRUE)) {
    ctl <- if(degree==0)list(regtype='lc') else if(degree==1)list(regtype='ll') else
      list(regtype='lp',degree=degree,bernstein.basis=bernstein)
    # Translated raw higher polynomials retain their established admission.
    shift <- if(bernstein)10000 else 0
    idx <- c(12,3,7,2,1,10,6,9,5,11,4,8)
    x <- data.frame(x=shift+4*d$x[idx]);y <- 1024+2*d$y[idx]
    for(tree in list(FALSE,TRUE,'auto')) {
      options(np.tree=tree,np.macMseries.accelerate=bernstein)
      b <- do.call(npcdensbw,c(list(xdat=x,ydat=y,bws=c(3,8),
        bwtype='generalized_nn',bwmethod='cv.ls',cxkertype='epanechnikov',
        cykertype='gaussian',bandwidth.compute=FALSE),ctl))
      # MPI_EVAL
      command <- substitute({
        options(np.tree=TREE,np.macMseries.accelerate=ACC,np.messages=FALSE)
        get('.npcdensbw_eval_only',asNamespace('npRmpi'))(X,Y,B,invalid.penalty='dbmax',force.local=FALSE)
      },list(TREE=tree,ACC=bernstein,X=x,Y=y,B=b))
      observed <- get('.npRmpi_bcast_cmd_expr',asNamespace('npRmpi'))(command,comm=1L,caller.execute=TRUE)$objective
      expect_lte(abs(observed-ref[[degree+1L]]['score']/2),1e-9)
      expect_true(is.finite(observed))
    }
  }
})

test_that('conditional GNN higher-order prefixes retain the literal criterion', {
  old <- options(np.messages=FALSE);on.exit(options(old),add=TRUE)
  skip_if_not(spawn_mpi_slaves(1L), "MPI session unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  d <- cgnn_fixture();x <- data.frame(x=d$x)
  expected <- cgnn_literal(d$x,d$y,8L,3L,2L,xorder=8L)['score']
  for(tree in c(FALSE,TRUE)) {
    options(np.tree=tree)
    b <- npcdensbw(xdat=x,ydat=d$y,bws=c(3,8),bwtype='generalized_nn',
      bwmethod='cv.ls',regtype='lp',degree=2L,bernstein.basis=TRUE,
      cxkertype='epanechnikov',cxkerorder=8L,cykertype='gaussian',
      bandwidth.compute=FALSE)
    command <- substitute({
      options(np.tree=TREE,np.messages=FALSE)
      get('.npcdensbw_eval_only',asNamespace('npRmpi'))(X,Y,B,
        invalid.penalty='dbmax',force.local=FALSE)
    },list(TREE=tree,X=x,Y=d$y,B=b))
    observed <- get('.npRmpi_bcast_cmd_expr',asNamespace('npRmpi'))(
      command,comm=1L,caller.execute=TRUE)$objective
    expect_lte(abs(observed-expected),1e-9)
  }
})

test_that('conditional GNN prefix tree admission cannot leak into other routes', {
  old <- options(np.messages=FALSE,np.tree=TRUE);on.exit(options(old),add=TRUE)
  skip_if_not(spawn_mpi_slaves(1L), "MPI session unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  d <- cgnn_fixture()
  b <- npcdensbw(xdat=data.frame(x=d$x),ydat=d$y,bws=c(3,8),bwtype='generalized_nn',
    bwmethod='cv.ls',cxkertype='epanechnikov',cykertype='gaussian',bandwidth.compute=FALSE)
  code <- get('.npcdensbw_tree_code',asNamespace('npRmpi'))
  yes <- get('DO_TREE_YES',asNamespace('npRmpi'))
  no <- get('DO_TREE_NO',asNamespace('npRmpi'))
  expect_identical(code(b,2L,0L),yes)
  uniform <- b;uniform$cxkertype <- 'uniform'
  expect_identical(code(uniform,2L,0L),yes)
  for(field in c('xncon','yncon','xnuno','xnord','ynuno','ynord')) {
    other <- b;other[[field]] <- other[[field]]+1L
    expect_identical(code(other,3L,1L),
                     if(field %in% c('xncon','xnuno','xnord'))yes else no)
  }
  for(order in c(4L,6L,8L)) {
    other <- b;other$cxkerorder <- order;expect_identical(code(other,2L,0L),yes)
  }
  for(order in c(4L,6L,8L)) {
    other <- b;other$cykerorder <- order
    expect_identical(code(other,2L,0L),yes)
    expect_identical(code(other,2L,0L,fit.context=TRUE),no)
  }
  other <- b;other$cykerorder <- 10L;expect_identical(code(other,2L,0L),no)
  other <- uniform;other$cxkerorder <- 4L;expect_identical(code(other,2L,0L),no)
  for(field in c('cxkerbound','cykerbound')) {
    other <- b;other[[field]] <- 'fixed';expect_identical(code(other,2L,0L),no)
  }
  other <- b;other$cxkertype <- 'gaussian';expect_identical(code(other,2L,0L),yes)
  for(kernel in c('epanechnikov','uniform')){
    other<-b;other$cykertype<-kernel
    expect_identical(code(other,2L,0L),yes)
  }
  other<-b;other$cykertype<-'uniform';other$cykerorder<-4L
  expect_identical(code(other,2L,0L),no)
  other<-b;other$cykertype<-'truncated gaussian'
  expect_identical(code(other,2L,0L),no)
  expect_identical(code(b,2L,0L,fit.context=TRUE),no)
  for(type in c('fixed','adaptive_nn')) {
    other <- b;other$type <- type;expect_identical(code(other,2L,0L),yes)
  }
  other <- b;other$method <- 'cv.ml';expect_identical(code(other,2L,0L),yes)
  options(np.tree=FALSE);expect_identical(code(b,2L,0L),no)
  expect_identical(code(uniform,2L,0L),no)
  options(np.tree='auto');expect_identical(code(b,2L,0L),no)
  expect_identical(code(uniform,2L,0L),no)
})

test_that('conditional GNN prefix preserves ties and extended-count geometry', {
  old<-options(np.messages=FALSE,np.extendednn=TRUE);on.exit(options(old),add=TRUE)
  skip_if_not(spawn_mpi_slaves(1L), "MPI session unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  ns<-asNamespace('npRmpi');ev<-get('.npcdensbw_eval_only',ns)
  for(tied in c(FALSE,TRUE)){
    d<-cgnn_fixture();if(tied){d$x[2]<-d$x[1];d$y[2]<-d$y[1]}
    k<-if(tied)c(3L,8L)else c(20L,20L)
    expected<-cgnn_literal(d$x,d$y,k[2],k[1],2L)['score']
    for(tree in c(FALSE,TRUE)){
      options(np.tree=tree)
      b<-npcdensbw(xdat=data.frame(x=d$x),ydat=d$y,bws=k,bwtype='generalized_nn',
        bwmethod='cv.ls',regtype='lp',degree=2L,bernstein.basis=TRUE,
        cxkertype='epanechnikov',cykertype='gaussian',bandwidth.compute=FALSE)
      observed<-ev(data.frame(x=d$x),d$y,b,invalid.penalty='dbmax')$objective
      expect_lte(abs(observed-expected),1e-9)
    }
  }
})

test_that('conditional GNN uniform prefixes preserve strict support and deletion', {
  old<-options(np.messages=FALSE,np.extendednn=TRUE);on.exit(options(old),add=TRUE)
  skip_if_not(spawn_mpi_slaves(1L), "MPI session unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  ns<-asNamespace('npRmpi')
  probe<-function(x,y,b,tree,acc){
    options(np.messages=FALSE,np.tree=tree,np.macMseries.accelerate=acc)
    get('.npcdensbw_eval_only',asNamespace('npRmpi'))(data.frame(x=x),y,b,
      invalid.penalty='dbmax',force.local=FALSE)$objective
  }
  for(kind in c('ordinary','ties','extended')){
    d<-cgnn_fixture()
    if(kind=='ties'){d$x[2]<-d$x[1];d$y[2]<-d$y[1]}
    k<-if(kind=='extended')c(20L,20L)else c(3L,8L)
    degrees<-if(kind=='ordinary')0:3 else 2L
    for(degree in degrees){
      expected<-cgnn_literal(d$x,d$y,k[2],k[1],degree,'uniform')['score']
      for(bernstein in c(FALSE,TRUE))for(tree in c(FALSE,TRUE)){
        options(np.tree=tree,np.macMseries.accelerate=bernstein)
        idx<-c(12,3,7,2,1,10,6,9,5,11,4,8)
        # Power-of-two scaling retains literal kth-distance ties.
        x<-4*d$x[idx];y<-2*d$y[idx]
        ctrl<-if(degree==0)list(regtype='lc')else if(degree==1)list(regtype='ll')else
          list(regtype='lp',degree=degree,bernstein.basis=bernstein)
        b<-do.call(npcdensbw,c(list(xdat=data.frame(x=x),ydat=y,bws=k,
          bwtype='generalized_nn',bwmethod='cv.ls',cxkertype='uniform',
          cykertype='gaussian',bandwidth.compute=FALSE),ctrl))
        cmd<-substitute(PROBE(X,Y,B,TREE,ACC),list(PROBE=probe,X=x,Y=y,B=b,TREE=tree,ACC=bernstein))
        observed<-get('.npRmpi_bcast_cmd_expr',ns)(cmd,comm=1L,caller.execute=TRUE)
        expect_true(is.finite(observed))
        expect_lte(abs(observed-expected/2),1e-9)
      }
    }
  }
})

test_that('conditional GNN prepared degree changes keep the canonical owner', {
  old<-options(np.messages=FALSE,np.tree=TRUE);on.exit(options(old),add=TRUE)
  skip_if_not(spawn_mpi_slaves(1L), "MPI session unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  d<-cgnn_fixture();x<-data.frame(x=d$x)
  ns<-asNamespace('npRmpi')
  b<-npcdensbw(xdat=x,ydat=d$y,bws=c(3,8),bwtype='generalized_nn',bwmethod='cv.ls',
    regtype='lp',degree=2L,bernstein.basis=TRUE,cxkertype='epanechnikov',
    cykertype='gaussian',bandwidth.compute=FALSE)
  probe<-function(x,y,b){
    ns<-asNamespace('npRmpi');options(np.messages=FALSE,np.tree=TRUE)
    prep<-get('.npcdensbw_prepared_prepare_args',ns)(xdat=x,ydat=y,bws=b,
      invalid.penalty='baseline',degree.search=TRUE)
    names(prep)[names(prep)=='penalty_mode']<-'penalty.mode'
    names(prep)[names(prep)=='penalty_multiplier']<-'penalty.multiplier'
    prepare<-get('npRmpiPreparedObjectivePrepareConditionalDensity',ns)
    evaluate<-get('npRmpiPreparedObjectiveEvalConditionalDensityRaw',ns)
    destroy<-get('npRmpiPreparedObjectiveDestroyConditionalDensity',ns)
    stopifnot(isTRUE(do.call(prepare,prep)))
    out<-tryCatch(vapply(0:3,function(q)evaluate(c(8,3),as.integer(q))[1],0),finally=destroy())
    stopifnot(isTRUE(do.call(prepare,prep)))
    again<-tryCatch(evaluate(c(8,3),2L)[1],finally=destroy())
    list(out=out,again=again)
  }
  for(kernel in c('epanechnikov','uniform')){
    b$cxkertype<-kernel
    # Reconstruct, so native kernel codes and public metadata agree.
    b<-npcdensbw(xdat=x,ydat=d$y,bws=c(3,8),bwtype='generalized_nn',bwmethod='cv.ls',
      regtype='lp',degree=2L,bernstein.basis=TRUE,cxkertype=kernel,
      cykertype='gaussian',bandwidth.compute=FALSE)
    expected<-vapply(0:3,function(q)cgnn_literal(d$x,d$y,8L,3L,q,kernel)['score'],0)
    command<-substitute(PROBE(X,Y,B),list(PROBE=probe,X=x,Y=d$y,B=b))
    observed<-get('.npRmpi_bcast_cmd_expr',ns)(command,comm=1L,caller.execute=TRUE)
    expect_equal(observed$out,expected,tolerance=1e-9)
    expect_equal(observed$again,expected[3],tolerance=1e-9)
  }
})

test_that('conditional GNN Gaussian response orders share the literal integral', {
  old<-options(np.messages=FALSE);on.exit(options(old),add=TRUE)
  skip_if_not(spawn_mpi_slaves(1L), "MPI session unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  d<-cgnn_fixture();x<-data.frame(x=d$x)
  for(yorder in c(4L,6L,8L))for(xorder in c(2L,4L,6L,8L,0L)){
    xkernel<-if(xorder==0L)'uniform'else'epanechnikov'
    xo<-if(xorder==0L)2L else xorder
    for(degree in c(0L,2L)){
      expected<-cgnn_literal(d$x,d$y,8L,3L,degree,xkernel,xo,yorder)['score']
      for(tree in c(FALSE,TRUE)){
        options(np.tree=tree)
        ctrl<-if(degree==0L)list(regtype='lc')else
          list(regtype='lp',degree=degree,bernstein.basis=TRUE)
        b<-do.call(npcdensbw,c(list(xdat=x,ydat=d$y,bws=c(3,8),
          bwtype='generalized_nn',bwmethod='cv.ls',cxkertype=xkernel,
          cykertype='gaussian',cykerorder=yorder,
          bandwidth.compute=FALSE),
          if(xorder==0L)list()else list(cxkerorder=xo),ctrl))
        command<-substitute({
          options(np.tree=TREE,np.messages=FALSE)
          get('.npcdensbw_eval_only',asNamespace('npRmpi'))(X,Y,B,
            invalid.penalty='dbmax',force.local=FALSE)
        },list(TREE=tree,X=x,Y=d$y,B=b))
        observed<-get('.npRmpi_bcast_cmd_expr',asNamespace('npRmpi'))(
          command,comm=1L,caller.execute=TRUE)$objective
        expect_lte(abs(observed-expected),1e-9)
      }
    }
  }
})

test_that('conditional GNN compact responses integrate the literal deleted criterion', {
  old<-options(np.messages=FALSE);on.exit(options(old),add=TRUE)
  skip_if_not(spawn_mpi_slaves(1L), "MPI session unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  d<-cgnn_fixture();x<-data.frame(x=d$x)
  for(yorder in c(2L,4L,6L,8L,0L))for(xorder in c(2L,4L,6L,8L,0L)){
    xkernel<-if(xorder==0L)'uniform'else'epanechnikov'
    ykernel<-if(yorder==0L)'uniform'else'epanechnikov'
    for(degree in c(0L,2L)){
      expected<-cgnn_literal(d$x,d$y,8L,3L,degree,xkernel,
        if(xorder==0L)2L else xorder,if(yorder==0L)2L else yorder,ykernel)['score']
      for(tree in c(FALSE,TRUE)){
        options(np.tree=tree)
        ctrl<-if(degree==0L)list(regtype='lc')else
          list(regtype='lp',degree=2L,bernstein.basis=TRUE)
        b<-do.call(npcdensbw,c(list(xdat=x,ydat=d$y,bws=c(3,8),
          bwtype='generalized_nn',bwmethod='cv.ls',cxkertype=xkernel,
          cykertype=ykernel,bandwidth.compute=FALSE),
          if(xorder==0L)list()else list(cxkerorder=xorder),
          if(yorder==0L)list()else list(cykerorder=yorder),ctrl))
        command<-substitute({
          options(np.tree=TREE,np.messages=FALSE)
          get('.npcdensbw_eval_only',asNamespace('npRmpi'))(X,Y,B,
            invalid.penalty='dbmax',force.local=FALSE)
        },list(TREE=tree,X=x,Y=d$y,B=b))
        observed<-get('.npRmpi_bcast_cmd_expr',asNamespace('npRmpi'))(
          command,comm=1L,caller.execute=TRUE)$objective
        expect_lte(abs(observed-expected),1e-9)
      }
    }
  }
})
