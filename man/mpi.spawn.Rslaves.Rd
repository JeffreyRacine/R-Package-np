\name{mpi.spawn.Rslaves}
\alias{mpi.spawn.Rslaves}
\alias{mpi.close.Rslaves}
\alias{tailslave.log}
\title{Spawn and Close R Slaves}
\description{
\code{mpi.spawn.Rslaves} spawns R slaves to those hosts automatically 
chosen by MPI or specific hosts assigned by the argument \code{hosts}.
Those R slaves are running in R BATCH mode with a specific Rscript file. 
The default Rscript file "slavedaemon.R" provides interactive R slave 
environments.

\code{mpi.close.Rslaves} shuts down R slaves spawned by 
\code{mpi.spawn.Rslaves}.

\code{tailslave.log} view (from tail) R slave log files (assuming they are all 
in one working directory).
}
\usage{
mpi.spawn.Rslaves(Rscript=system.file("slavedaemon.R", package="npRmpi"),
                  nslaves=mpi.universe.size(),
                  root = 0,
                  intercomm = 2,
                  comm = 1,
                  hosts = NULL,
                  needlog = FALSE,
                  mapdrive=TRUE,
                  quiet = FALSE,
                  nonblock=TRUE,
                  sleep=0.1)

mpi.close.Rslaves(dellog = TRUE, comm = 1, force = FALSE)
tailslave.log(nlines = 3, comm = 1)
}
\arguments{
  \item{Rscript}{an R script file used to run R in BATCH mode.}
  \item{nslaves}{number of slaves to be spawned.}
  \item{root}{the rank number of the member who spawns R slaves.}
  \item{intercomm}{an intercommunicator number}
  \item{comm}{a communicator number merged from an intercomm.}
  \item{hosts}{NULL or LAM node numbers to specify where R slaves are to be 
spawned.}
  \item{needlog}{a logical. If TRUE, R BATCH outputs will be saved in log files.  
  If FALSE, the outputs will send to /dev/null.}
  \item{mapdrive}{a logical. If TRUE and master's working dir is on a network, mapping network drive 
is attemped on remote nodes under windows platform.}
  \item{quiet}{a logical. If TRUE, do not print anything unless an error occurs.}
  \item{nonblock}{a logical. If TRUE, a nonblock procedure is used on all slaves so that
they will consume none or little CPUs while waiting.}
  \item{sleep}{a sleep interval, used when nonblock=TRUE. The smaller sleep is, the more responsive slaves
are, the more CPUs consume.}
  \item{dellog}{a logical specifying if R slave's log files are deleted or 
not.}
  \item{force}{a logical. If \code{TRUE}, force a hard shutdown of slave daemons.
  When \code{options(npRmpi.reuse.slaves=TRUE)} and \code{force=FALSE},
  \code{mpi.close.Rslaves()} performs a soft-close (i.e., keeps daemons alive for reuse).}
  \item{nlines}{number of lines to view from tail in R slave's log files.}
}
\value{
 Unless \code{quiet = TRUE}, \code{mpi.spawn.Rslaves} prints to stdio how many 
slaves are successfully spawned and where they are running.

 \code{mpi.close.Rslaves} returns a status code. When \code{options(npRmpi.reuse.slaves=TRUE)}
 and \code{force=FALSE}, this may be a no-op (soft-close) so that spawned daemons
 can be reused within the same R session.

\code{tailslave.log} returns last lines of R slave's log files.
}
\details{
The R slaves that \code{mpi.spawn.Rslaves} spawns are really running a shell 
program which can be found in \code{system.file("Rslaves.sh",package="npRmpi")} 
which takes a Rscript file as one of its arguments. Other arguments are 
used to see if a log file (R output) is needed and how to name it. The master 
process id and the comm number, along with host names where R slaves are running  are 
used to name these log files. 

Once R slaves are successfully spawned, the 
mergers from an intercomm (default `intercomm = 2') to a comm (default 
`comm = 1') are automatically done on master and slaves (should be done 
if the default Rscript is replaced). If additional sets of R slaves are 
needed, please use `comm = 3', `comm = 4', etc to spawn them. At most a 
comm number up to 10 can be used. Notice that the default comm number for 
R slaves (using slavedaemon.R) is always 1 which is saved as .comm.  

On some systems (notably macOS+MPICH), repeatedly spawning and tearing down
slaves in the same R session can lead to hangs/crashes. To avoid this,
npRmpi may reuse an existing slave pool when \code{options(npRmpi.reuse.slaves=TRUE)}.
In this mode, \code{mpi.spawn.Rslaves()} becomes idempotent and
\code{mpi.close.Rslaves(force=FALSE)} performs a soft-close.
  
To spawn R slaves to specific hosts, please use the argument \code{hosts} 
with a list of those node numbers (an integer vector). Total node numbers 
along their host names can be found by using \code{\link{mpi.hostinfo}}.
Notice that this is MPI implementation specific. 
} 

\seealso{
  \code{\link{mpi.comm.spawn}},
  \code{\link{mpi.hostinfo}}.
}
\author{
Hao Yu
}
\examples{
\dontrun{
# Not run in checks: spawning/tearing down MPI daemons is environment-dependent
# and can interfere with later examples in the same session.
mpi.spawn.Rslaves(nslaves=2)
tailslave.log()
mpi.remote.exec(rnorm(10))
mpi.close.Rslaves()
}
}
\keyword{utilities}
