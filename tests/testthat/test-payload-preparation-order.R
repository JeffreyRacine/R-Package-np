test_that('payload preparation errors do not schedule worker receives', {
ns<-asNamespace('npRmpi');raw<-list()
roles<-c('invocation','object','all_objects','all_functions','scatter','remote','parsim')
selected<-c('.npRmpi_invocation_bcast','mpi.bcast.Robj','mpi.bcast.Robj2slave','mpi.bcast.Rfun2slave',
 '.npRmpi_bcast_object_to_slaves','mpi.scatter.Robj2slave','mpi.scatter.Robj','.npRmpi_scatter_prepare',
 'mpi.remote.exec','mpi.parSim')
for(role in roles){
 e<-new.env(parent=ns)
 for(name in selected)if(exists(name,ns,inherits=FALSE)) {
  fun<-get(name,ns);environment(fun)<-e;e[[name]]<-fun
 }
 events<-character();frame<-new.env(parent=emptyenv())
 err<-structure(list(message='INJECTED serialization failure',call=quote(serialize(packet,NULL)),code=171L),class=c('preparation_error','error','condition'))
 e$serialize<-function(object,connection,...){events<<-c(events,'serialize');stop(err)}
 e$mpi.comm.rank<-function(...)0L;e$mpi.comm.size<-function(...)2L
 e$mpi.bcast.cmd<-function(...){events<<-c(events,'schedule');invisible(NULL)}
 e$.npRmpi_bcast_cmd_expr<-e$mpi.bcast.cmd
 e$.npRmpi_invocation_frame<-function(...)frame
 e$.npRmpi_bcast_prepared<-function(...){events<<-c(events,'transfer');invisible(NULL)}
 e$.npRmpi_transport_trace<-function(...)invisible(NULL)
 e$objects<-function(...)'bounded_payload'
 e$mget<-function(...)list(bounded_payload=if(role=='all_functions')function(x)x+1 else 42L)
 e$payload<-42L
 call<-switch(role,
  invocation=quote(.npRmpi_invocation_bcast('test','payload',payload)),
  object=quote(mpi.bcast.Robj2slave(payload)),
  all_objects=quote(mpi.bcast.Robj2slave(all=TRUE)),
  all_functions=quote(mpi.bcast.Rfun2slave()),
  scatter=quote(mpi.scatter.Robj2slave(list(payload))),
  remote=quote(mpi.remote.exec(sum,1:3)),
  parsim=quote(mpi.parSim(n=3,nsim=1,statistic=mean)))
 e$runif<-function(...)0.5
 got<-tryCatch(eval(call,e),error=identity)
 okay<-identical(got,err)&&identical(events,'serialize')&&!length(ls(frame,all.names=TRUE))
 raw[[role]]<-list(condition=got,events=events,frame=ls(frame,all.names=TRUE),okay=okay)
 expect_true(okay,info=paste(role,paste(events,collapse=' -> ')))
}
})
