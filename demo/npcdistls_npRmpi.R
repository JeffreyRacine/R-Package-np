## Make sure you have the .Rprofile file from npRmpi/inst/ in your
## current directory or home directory. It is necessary.

## To run this on systems with OPENMPI installed and working, try
## mpirun -np 2 R CMD BATCH npcdistml_npRmpi. Check the time in the
## output file foo.Rout (the name of this file with extension .Rout),
## then try with, say, 4 processors and compare run time.

## Initialize master and slaves.

npRmpi.start(nslaves=1)

## Turn off progress i/o as this clutters the output file (if you want
## to see search progress you can comment out this command)

options(npRmpi.autodispatch=TRUE, np.messages=FALSE)

## Load your data and broadcast it to all slave nodes

library(MASS)

set.seed(42)

n <- 2000

rho <- 0.25
mu <- c(0,0)
Sigma <- matrix(c(1,rho,rho,1),2,2)
data <- mvrnorm(n=n, mu, Sigma)
mydat <- data.frame(x=data[,2],y=data[,1])
## A conditional distribution estimation example. 

t <- system.time(bw <- npcdistbw(y~x,
                                               bwmethod="cv.ls",
                                               data=mydat))

summary(bw)

t <- t + system.time(model <- npcdist(bws=bw))

summary(model)

cat("Elapsed time =", t[3], "\n")

## Clean up properly then quit()

npRmpi.stop(force=TRUE)
