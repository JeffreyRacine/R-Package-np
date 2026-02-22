library(npRmpi)

## To run this on systems with OPENMPI installed and working, try
## mpirun -np 2 R CMD BATCH npcmstest_npRmpi. Check the time in the
## output file foo.Rout (the name of this file with extension .Rout),
## then try with, say, 4 processors and compare run time.

## Initialize master and slaves.

## Batch/cluster usage (attach mode under mpiexec):
##   mpiexec -n <master+slaves> R CMD BATCH --vanilla <script>.R
## Inside the script, use attach mode instead of spawning:
##   try(mpi.comm.dup(0, 1), silent = TRUE)
##   npRmpi.init(mode="attach", comm=1, autodispatch=TRUE, np.messages=FALSE)
##
npRmpi.init(mode="attach", comm=1, autodispatch=TRUE)

## Turn off progress i/o as this clutters the output file (if you want
## to see search progress you can comment out this command)
data("oecdpanel")
attach(oecdpanel)

n <- 616 # n <- nrow(oecdpanel)

oecd <- factor(oecd)
year <- factor(year)

model <- lm(growth ~ oecd +
            year +
            initgdp +
            I(initgdp^2) +
            I(initgdp^3) +
            I(initgdp^4) +
            popgro +
            inv +
            humancap +
            I(humancap^2) +
            I(humancap^3) - 1, 
            x=TRUE, 
            y=TRUE)

X <- data.frame(oecd, year, initgdp, popgro, inv, humancap)
## Consistent model specification test (we override defaults for
## demonstration purposes - don't do this for real problems).

t <- system.time(output <- npcmstest(model = model,
                                                   xdat = X,
                                                   ydat = growth,
                                                   nmulti=1,
                                                   ftol=.01,
                                                   tol=.01))

output

cat("Elapsed time =", t[3], "\n")

## Clean up properly then quit()

npRmpi.quit(mode="attach", comm=1)
## Batch/cluster attach-mode shutdown (for mpiexec workflows):
##   npRmpi.quit(mode="attach", comm=1)
## (no force=TRUE required for attach mode)
