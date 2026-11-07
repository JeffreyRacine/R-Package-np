test_that('required retained training roles cannot recurse through NULL', {
  withr::local_options(np.messages=FALSE)
  ns<-asNamespace('npRmpi');resolve<-get('.np_retained_training_args',ns)
  specs<-readRDS(test_path('fixtures','required-training-specs.rds'));required.count<-0L;optional.count<-0L
  for(p in specs) {
    definition<-get(p[[1]],ns);roles<-p[[2]];training<-if(length(p)>2)p[[3]] else names(roles)
    args<-setNames(rep(list(1:5),length(roles)),names(roles))
    expect_identical(resolve(list(),roles,args,definition,training),c(list(bws=list()),args))
    for(arg in training) {
      a<-args;a[arg]<-list(NULL);default<-formals(definition)[[arg]]
      if(is.call(default)&&identical(default[[1L]],quote(stop))) {
        required.count<-required.count+1L
        expect_error(resolve(list(),roles,a,definition,training),paste0("training data '",arg,"' missing"),fixed=TRUE)
      } else {
        optional.count<-optional.count+1L
        expect_identical(resolve(list(),roles,a,definition,training),c(list(bws=list()),a))
      }
    }
  }
  expect_identical(required.count,20L);expect_identical(optional.count,1L)
  count<-0L;null<-function(){count<<-count+1L;NULL}
  expect_error(resolve(list(),c(tdat='dat'),list(tdat=null()),get('npudens.bandwidth',ns)),"training data 'tdat' missing",fixed=TRUE)
  expect_identical(count,1L)
  # A NULL optional response on the hat route is not a required training role.
  args<-list(txdat=1:5,y=NULL)
  expect_identical(resolve(list(),c(txdat='xdat',y='ydat'),args,get('npreghat.rbandwidth',ns),training='txdat'),c(list(bws=list()),args))
})

test_that('public required-NULL calls fail promptly and valid calls still work', {
  if (exists('spawn_mpi_slaves', mode='function')) {
    spawn_mpi_slaves()
    on.exit(close_mpi_slaves(force=TRUE), add=TRUE)
  }
  withr::local_options(np.messages=FALSE)
  n<-24L;x<-data.frame(x=seq(.03,.97,length.out=n));z<-data.frame(z=sin(seq_len(n)));y<-x$x^2+.03*cos(seq_len(n));yd<-data.frame(y=y)
  rb<-npregbw(xdat=x,ydat=y,bws=.2,bandwidth.compute=FALSE,regtype='ll')
  cb<-npcdensbw(xdat=x,ydat=yd,bws=c(.2,.2),bandwidth.compute=FALSE)
  db<-npcdistbw(xdat=x,ydat=yd,bws=c(.2,.2),bandwidth.compute=FALSE)
  mb<-npcdensbw(xdat=x,ydat=data.frame(g=factor(rep(c('a','b'),n/2))),bws=c(.2,.1),bandwidth.compute=FALSE)
  bw<-list(npudens=npudensbw(dat=x,bws=.2,bandwidth.compute=FALSE),npudist=npudistbw(dat=x,bws=.2,bandwidth.compute=FALSE),npreg=rb,npcdens=cb,npcdist=db,npindex=npindexbw(xdat=cbind(x,z),ydat=y,bws=c(1,.3,.3),method='ichimura',bandwidth.compute=FALSE),npscoef=npscoefbw(xdat=x,ydat=y,zdat=z,bws=.2,bandwidth.compute=FALSE),npplreg=npplregbw(xdat=x,ydat=y,zdat=z,bws=matrix(c(.2,.2),2),bandwidth.compute=FALSE),npqreg=db,npconmode=mb,npreghat=rb)
  specs<-readRDS(test_path('fixtures','required-training-specs.rds'));ns<-asNamespace('npRmpi')
  for(family in names(specs)) {
    p<-specs[[family]];training<-if(length(p)>2)p[[3]] else names(p[[2]]);definition<-get(p[[1]],ns)
    for(arg in training) {
      default<-formals(definition)[[arg]]
      if(is.call(default)&&identical(default[[1L]],quote(stop))) {
        a<-c(list(bws=bw[[family]]),setNames(list(NULL),arg))
        expect_error(do.call(get(family),a),paste0("training data '",arg,"' missing"),fixed=TRUE)
      }
    }
  }
  # Required failure leaves ordinary fitting usable, including response-only refit.
  a<-npreg(rb,se=FALSE);b<-npreg(rb,tydat=y+1,se=FALSE)
  expect_equal(fitted(b),fitted(a)+1,tolerance=1e-12)
  expect_true(all(is.finite(fitted(npudens(bw$npudens)))))
  expect_true(all(is.finite(fitted(npudist(bw$npudist)))))
})
