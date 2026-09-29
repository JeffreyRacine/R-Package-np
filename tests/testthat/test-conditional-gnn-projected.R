# Independent deleted-WLS/full-line oracle; no native radius or row helper.
cgnn_projected_fixture <- function() list(
  x=c(-.67932861763984,.597772665787488,-.135870280675590,-.180028952658176,
      -.0812444845214486,-.0922158942557871,-.119163474533707,.565201563760638,
      -.167390177026391,-.617473398335278,-.452014316339046,.213771468494087),
  y=c(.117393371024781,-.215783902381389,-.752371990755028,1.55294725882662,
      .776713307792541,.329678953352629,-1.89405512568533,-1.99349412668255,
      -.220072074888336,.198501257716297,.563714301921436,-.880794421005045))

cgnn_projected_literal<-function(d,degree,kx,ky){
 n<-length(d$x);B<-outer(d$x,0:degree,`^`)
 radius<-function(a,z,k){i<-min(k,length(a)-1L);sort(abs(a-z))[i]*k/i}
 rows<-lapply(seq_len(n),function(i){
  w<-dnorm((d$x[i]-d$x[-i])/radius(d$x[-i],d$x[i],kx));A<-B[-i,,drop=FALSE]
  w*as.vector(A%*%solve(crossprod(A,w*A),B[i,]))
 })
 cross<-vapply(seq_len(n),function(i){h<-radius(d$y[-i],d$y[i],ky);sum(rows[[i]]*dnorm((d$y[i]-d$y[-i])/h)/h)},0)
 cuts<-c(-Inf,sort(unique(c(d$y,as.vector(outer(d$y,d$y,`+`)/2)))),Inf)
 I1<-sum(vapply(seq_len(length(cuts)-1L),function(j)
  integrate(function(q)vapply(q,function(z)mean(vapply(seq_len(n),function(i){h<-radius(d$y[-i],z,ky);sum(rows[[i]]*dnorm((z-d$y[-i])/h)/h)^2},0)),0),
    cuts[j],cuts[j+1L],abs.tol=1e-12/length(cuts),rel.tol=1e-10,subdivisions=1000L)$value,0))
 list(I1=I1,I2=mean(cross),cross=cross,score=2*mean(cross)-I1)
}

test_that('conditional Gaussian GNN projected rows retain the literal criterion', {
  old<-options(np.messages=FALSE);on.exit(options(old),add=TRUE)
  skip_if_not(spawn_mpi_slaves(1L), 'MPI session unavailable')
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  d<-cgnn_projected_fixture();x<-data.frame(x=d$x);y<-d$y
  for(degree in 0:3){
    ref<-cgnn_projected_literal(d,degree,8L,3L)
    for(tree in c(FALSE,TRUE)){
      options(np.tree=tree)
      ctl<-if(degree==0)list(regtype='lc')else if(degree==1)list(regtype='ll')else
        list(regtype='lp',degree=degree,bernstein.basis=TRUE)
      b<-do.call(npcdensbw,c(list(xdat=x,ydat=y,bws=c(3,8),bwtype='generalized_nn',
        bwmethod='cv.ls',bandwidth.compute=FALSE),ctl))
      observed<-{
    command<-substitute({options(np.messages=FALSE,np.tree=TREE);
      get('.npcdensbw_eval_only',asNamespace('npRmpi'))(X,Y,B,invalid.penalty='dbmax',force.local=FALSE)$objective},
      list(TREE=tree,X=x,Y=y,B=b))
    get('.npRmpi_bcast_cmd_expr',asNamespace('npRmpi'))(command,comm=1L,caller.execute=TRUE)
  }
      expect_lte(abs(observed-ref$score),1e-9)
      expect_true(is.finite(observed))
    }
  }
})
