## This is the serial version of npsigtest_npRmpi.R for comparison
## purposes (bandwidth ought to be identical, timing may
## differ). Study the differences between this file and its MPI
## counterpart for insight about your own problems.

library(np)
options(np.messages=FALSE)

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
