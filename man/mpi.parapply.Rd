\name{mpi.applyLB}
\alias{mpi.applyLB}
\alias{mpi.parApply}
\alias{mpi.parLapply}
\alias{mpi.parSapply}
\alias{mpi.parRapply}
\alias{mpi.parCapply}
\alias{mpi.parReplicate}
\alias{mpi.parMM}
\title{(Load balancing) parallel apply}
\description{
(Load balancing) parallel\code{lapply} and related functions. 
}
\usage{
mpi.applyLB(X, FUN, ..., apply.seq=NULL, comm=1)
mpi.parApply(X, MARGIN, FUN, ..., job.num = mpi.comm.size(comm)-1,
                    apply.seq=NULL, comm=1)
mpi.parLapply(X, FUN, ..., job.num=mpi.comm.size(comm)-1, apply.seq=NULL, 
		comm=1)  
mpi.parSapply(X, FUN, ..., job.num=mpi.comm.size(comm)-1, apply.seq=NULL, 
		simplify=TRUE, USE.NAMES = TRUE, comm=1)  
mpi.parRapply(X, FUN, ..., job.num=mpi.comm.size(comm)-1, apply.seq=NULL, 
		comm=1)  
mpi.parCapply(X, FUN, ..., job.num=mpi.comm.size(comm)-1, apply.seq=NULL, 
		comm=1)  
mpi.parReplicate(n, expr, job.num=mpi.comm.size(comm)-1, apply.seq=NULL, 
		simplify = TRUE, comm=1)
mpi.parMM (A, B, job.num=mpi.comm.size(comm)-1, comm=1)
}
\arguments{
  \item{X}{an array or matrix.}
  \item{MARGIN}{vector specifying the dimensions to use.}
  \item{FUN}{a function.}
  \item{simplify}{logical; should the result be simplified to a vector or
          matrix if possible?}
  \item{USE.NAMES}{logical; if \code{TRUE} and if \code{X} is character, use \code{X} as
          \code{names} for the result unless it had names already.}
  \item{n}{number of replications.}
  \item{A}{a matrix}
  \item{B}{a matrix}
  \item{expr}{expression to evaluate repeatedly.}
  \item{job.num}{Total job numbers. If job numbers is bigger than total slave numbers (default 
value), a load balancing approach is used.}
  \item{apply.seq}{if reproducing the same computation (simulation) is desirable, set it 
to the integer vector .mpi.applyLB generated in previous computation (simulation).}
  \item{...}{optional arguments to \code{FUN}}
  \item{comm}{a communicator number}
}
\details{
Unless length of \code{X} is no more than total slave numbers (slave.num) and in this case 
\code{mpi.applyLB} is the same as \code{mpi.apply}, \code{mpi.applyLB} sends a next job to a 
slave who just delivered a finished job. The sequence of slaves who deliver results to master are 
saved into \code{.mpi.applyLB}. It keeps track which part of results done by which slaves. 
\code{.mpi.applyLB} can be used to reproduce the same simulation result if the same seed is 
used and the argument \code{apply.seq} is equal to \code{.mpi.applyLB}.
 
With the default value of argument \code{job.num} which is slave.num, \code{mpi.parApply}, 
\code{mpi.parLapply}, \code{mpi.parSapply}, \code{mpi.parRapply}, \code{mpi.parCapply}, 
\code{mpi.parSapply}, and \code{mpi.parMM} are clones of \pkg{snow}'s parApply, parLappy, 
parSapply, parRapply, parCapply, parSapply, and parMM, respectively. When \code{job.num} is 
bigger than slave.num, a load balancing approach is used.
 }
\section{Warning}{
When using the argument \code{apply.seq} with \code{.mpi.applyLB}, be sure all settings are the same 
as before, i.e., the same data, job.num, slave.num, and seed. Otherwise a deadlock could occur. 
Notice that \code{apply.seq} is useful only if \code{job.num} is bigger than slave.num. 
}
\seealso{
\code{\link{mpi.apply}}
}
\examples{
\donttest{
#Assume that there are some slaves running

#mpi.applyLB
x=1:7
mpi.applyLB(x,rnorm,mean=2,sd=4)

#get the same simulation 
mpi.remote.exec(set.seed(111))
mpi.applyLB(x,rnorm,mean=2,sd=4)
mpi.remote.exec(set.seed(111))
mpi.applyLB(x,rnorm,mean=2,sd=4,apply.seq=.mpi.applyLB)

#mpi.parApply
x=1:24
dim(x)=c(2,3,4)
mpi.parApply(x, MARGIN=c(1,2), FUN=mean,job.num = 5)

#mpi.parLapply
mdat <- matrix(c(1,2,3, 7,8,9), nrow = 2, ncol=3, byrow=TRUE,
                    dimnames = list(c("R.1", "R.2"), c("C.1", "C.2", "C.3")))
mpi.parLapply(mdat, rnorm) 

#mpi.parSapply
mpi.parSapply(mdat, rnorm) 

#mpi.parMM
A=matrix(1:1000^2,ncol=1000)
mpi.parMM(A,A)
}
}
\keyword{utilities}
