library(npRmpi)

## To run this on systems with OPENMPI installed and working, try
## mpirun -np 2 R CMD BATCH npconmode_npRmpi. Check the time in the
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
library(MASS)

data(birthwt)              

n <- 189 # n <- nrow(birthwt)

birthwt$low <- factor(birthwt$low)
birthwt$smoke <- factor(birthwt$smoke)
birthwt$race <- factor(birthwt$race)
birthwt$ht <- factor(birthwt$ht)
birthwt$ui <- factor(birthwt$ui)
birthwt$ftv <- ordered(birthwt$ftv)
## A conditional mode example

t <- system.time(bw <- npcdensbw(low~
                                               smoke+ 
                                               race+ 
                                               ht+ 
                                               ui+    
                                               ftv+  
                                               age+           
                                               lwt,
                                               data=birthwt))

summary(bw)

t <- t + system.time(model <- npconmode(bws=bw))

summary(model)

cat("Elapsed time =", t[3], "\n")

## Clean up properly then quit()

npRmpi.quit(mode="attach", comm=1)
## Batch/cluster attach-mode shutdown (for mpiexec workflows):
##   npRmpi.quit(mode="attach", comm=1)
## (no force=TRUE required for attach mode)
