## Profile/manual-broadcast demo (mpiexec + .Rprofile + mpi.bcast.*).
## Run with two ranks (master + one worker), e.g.
##   mpiexec -env R_PROFILE_USER ../.Rprofile -env R_PROFILE "" \\
##           -n 2 R CMD BATCH --no-save <script>.R
## Do not use R CMD BATCH --vanilla for profile mode.
##
## Initialize master and slaves.

mpi.bcast.cmd(np.mpi.initialize(),
              caller.execute=TRUE)

## Turn off progress i/o as this clutters the output file (if you want
## to see search progress you can comment out this command)

mpi.bcast.cmd(options(np.messages=FALSE),
              caller.execute=TRUE)

mpi.bcast.cmd(data("oecdpanel"),
              caller.execute=TRUE)
mpi.bcast.cmd(attach(oecdpanel),
              caller.execute=TRUE)

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

mpi.bcast.Robj2slave(model)
mpi.bcast.Robj2slave(X)

## Consistent model specification test (we override defaults for
## demonstration purposes - don't do this for real problems).

t <- system.time(mpi.bcast.cmd(output <- npcmstest(model = model,
                                                   xdat = X,
                                                   ydat = growth,
                                                   nmulti=1,
                                                   ftol=.01,
                                                   tol=.01),
                               caller.execute=TRUE))

output

cat("Elapsed time =", t[3], "\n")

## Clean up properly then quit()

mpi.bcast.cmd(mpi.quit(),
              caller.execute=TRUE)
