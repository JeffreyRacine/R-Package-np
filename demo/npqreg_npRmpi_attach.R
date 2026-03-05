library(npRmpi)

## Attach mode demo (session routing under mpiexec).
## Run with two ranks (master + one worker), e.g.
##   mpiexec -n 2 Rscript --no-save <script>.R
## or
##   mpiexec -n 2 R CMD BATCH --no-save <script>.R
##
## If running under mpiexec, keep profile env vars cleared for attach mode:
##   -env R_PROFILE_USER "" -env R_PROFILE ""
##
## Initialize master and slaves.
npRmpi.init(mode="attach", comm=1, autodispatch=TRUE)
options(np.messages=FALSE)

if (mpi.comm.rank(0L) == 0L) {

## Turn off progress i/o as this clutters the output file (if you want
## to see search progress you can comment out this command)
data("Italy")



n <- 1008 # n <- nrow(Italy)
## A quantile regression example

t <- system.time(bw <- npcdistbw(gdp~ordered(year),data=Italy))

summary(bw)

t <- t + system.time(model.q0.25 <- npqreg(bws=bw, tau=0.25))
t <- t + system.time(model.q0.50 <- npqreg(bws=bw, tau=0.50))
t <- t + system.time(model.q0.75 <- npqreg(bws=bw, tau=0.75))

cat("Elapsed time =", t[3], "\n")

## Clean up properly then quit()

npRmpi.quit(mode="attach", comm=1)
}
## Batch/cluster attach-mode shutdown (for mpiexec workflows):
##   npRmpi.quit(mode="attach", comm=1)
## (no force=TRUE required for attach mode)
