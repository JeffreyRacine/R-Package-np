## Make sure you have the .Rprofile file from npRmpi/inst/ in your
## current directory or home directory. It is necessary.

## To run this on systems with OPENMPI installed and working, try
## mpirun -np 2 R CMD BATCH npsigtest_npRmpi. Check the time in the
## output file foo.Rout (the name of this file with extension .Rout),
## then try with, say, 4 processors and compare run time.

## Initialize master and slaves.
npRmpi.start(nslaves=1)

## Turn off progress i/o as this clutters the output file (if you want
## to see search progress you can comment out this command)
options(npRmpi.autodispatch=TRUE, np.messages=FALSE)
set.seed(42)

## Significance testing with z irrelevant

n <- 1000

z <- factor(rbinom(n,1,.5))
x1 <- rnorm(n)
x2 <- runif(n,-2,2)
y <- x1 + x2 + rnorm(n)
mydat <- data.frame(z,x1,x2,y)
rm(z,x1,x2,y)

t <- system.time(model <- npreg(y~z+x1+x2,
                                regtype="ll",
                                bwmethod="cv.aic",
                                data=mydat))

## An example of the consistent nonparametric significance test

t <- t + system.time(output <- npsigtest(model))

output

cat("Elapsed time =", t[3], "\n")

## Clean up properly
npRmpi.stop(force=TRUE)
