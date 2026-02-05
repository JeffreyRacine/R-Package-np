\name{mpi.wait}
\title{Nonblocking completion operations}
\alias{mpi.cancel}
\alias{mpi.test.cancelled}
\alias{mpi.test}
\alias{mpi.testall}
\alias{mpi.testany}
\alias{mpi.testsome}
\alias{mpi.wait}
\alias{mpi.waitall}
\alias{mpi.waitany}
\alias{mpi.waitsome}
\description{
\code{mpi.cancel} cancels a nonblocking send or receive request.

\code{mpi.test.cancelled} tests if \code{mpi.cancel} cancels or not.

\code{wait}, \code{waitall}, \code{waitany}, and \code{waitsome} are used 
to complete nonblocking send or receive requests. They are not local.

\code{test}, \code{testall}, \code{testany}, and \code{testsome} are used 
to complete nonblocking send and receive requests. They are local.
}
\usage{
mpi.cancel(request)
mpi.test.cancelled(status=0)
mpi.test(request, status=0)
mpi.testall(count)
mpi.testany(count, status=0)
mpi.testsome(count)
mpi.wait(request, status=0)
mpi.waitall(count)
mpi.waitany(count, status=0)
mpi.waitsome(count)
}
\arguments{
  \item{count}{total number of nonblocking operations.}
  \item{request}{a request number.}
  \item{status}{a status number.}
}
\value{
\code{mpi.cancel} returns no value.

\code{mpi.test.cancelled} returns TRUE if a nonblocking call is cancelled; 
FALSE otherwise.

\code{mpi.wait} returns no value. Instead status contains information that 
can be retrieved by \code{mpi.get.count} and \code{mpi.get.sourcetag}.

\code{mpi.test} returns TRUE if a request is complete; FALSE otherwise. If 
TRUE, it is the same as \code{mpi.wait}.

\code{mpi.waitany} returns which request (index) has been completed. In
addition, status contains information that can be retrieved by
\code{mpi.get.count} and \code{mpi.get.sourcetag}. 

\code{mpi.testany} returns a list: index--- request index; flag---TRUE if 
a request is complete; FALSE otherwise (index is no use in this case). 
If flag is TRUE, it is the same as \code{mpi.waitany}.

\code{mpi.waitall} returns no value. Instead statuses 0, 1, ..., count-1 
contain corresponding information that can be retrieved by 
\code{mpi.get.count} and \code{mpi.get.sourcetag}.

\code{mpi.testall} returns TRUE if all requests are complete; FALSE 
otherwise. If TRUE, it is the same as \code{mpi.waitall}.

\code{mpi.waitsome} returns a list: count--- number of requests that have
been completed; indices---an integer vector of size \$count of those 
completed request numbers (in  0, 1 ,..., count-1). In addition, statuses 
0, 1, ..., \$count-1 contain corresponding information that can be 
retrieved by \code{mpi.get.count} and \code{mpi.get.sourcetag}.  

\code{mpi.testsome} is the same as \code{mpi.waitsome} except that \$count 
may be 0 and in this case \$indices is no use.
}
\details{
\code{mpi.wait} and \code{mpi.test} are used to complete a nonblocking 
send and receive request: use the same request number by \code{mpi.isend} 
or \code{mpi.irecv}. Once completed, the associated request is set to 
MPI\_REQUEST\_NULL and status contains information such as source, tag, 
and length of message.

If multiple nonblocking sends or receives are initiated, the following 
calls are more efficient. Make sure that request numbers are used 
consecutively as request=0, request=1, request=2, etc. In this way, the 
following calls can find request information in system memory.

\code{mpi.waitany} and \code{mpi.testany} are used to complete one out of 
several requests. 

\code{mpi.waitall} and \code{mpi.testall} are used to complete all 
requests. 

\code{mpi.waitsome} and \code{mpi.testsome} are used to complete all 
enabled requests. 
}
\author{
  Hao Yu
}
\seealso{
\code{\link{mpi.isend}},
\code{\link{mpi.irecv}},
\code{\link{mpi.get.count}},
\code{\link{mpi.get.sourcetag}}.
}
\references{
  \url{http://www.lam-mpi.org/},  
  \url{http://www-unix.mcs.anl.gov/mpi/www/www3/}
}
%\examples{
%}
\keyword{interface}
