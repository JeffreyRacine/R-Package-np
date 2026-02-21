## To run this on systems with OPENMPI installed and working, try
## mpirun -np 2 R CMD BATCH npconmode_npRmpi. Check the time in the
## output file foo.Rout (the name of this file with extension .Rout),
## then try with, say, 4 processors and compare run time.

## Initialize master and slaves.

npRmpi.start(nslaves=1)

## Turn off progress i/o as this clutters the output file (if you want
## to see search progress you can comment out this command)

options(npRmpi.autodispatch=TRUE, np.messages=FALSE)

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

npRmpi.stop(force=TRUE)
