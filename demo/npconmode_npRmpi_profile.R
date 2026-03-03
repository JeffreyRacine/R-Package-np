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

mpi.bcast.cmd(library(MASS),
              caller.execute=TRUE)

mpi.bcast.cmd(data(birthwt),
              caller.execute=TRUE)              

n <- 189 # n <- nrow(birthwt)

birthwt$low <- factor(birthwt$low)
birthwt$smoke <- factor(birthwt$smoke)
birthwt$race <- factor(birthwt$race)
birthwt$ht <- factor(birthwt$ht)
birthwt$ui <- factor(birthwt$ui)
birthwt$ftv <- ordered(birthwt$ftv)

mpi.bcast.Robj2slave(birthwt)

## A conditional mode example

t <- system.time(mpi.bcast.cmd(bw <- npcdensbw(low~
                                               smoke+ 
                                               race+ 
                                               ht+ 
                                               ui+    
                                               ftv+  
                                               age+           
                                               lwt,
                                               data=birthwt),
                               caller.execute=TRUE))

summary(bw)

t <- t + system.time(mpi.bcast.cmd(model <- npconmode(bws=bw),
                                   caller.execute=TRUE))

summary(model)

cat("Elapsed time =", t[3], "\n")

## Clean up properly then quit()

mpi.bcast.cmd(mpi.quit(),
              caller.execute=TRUE)
