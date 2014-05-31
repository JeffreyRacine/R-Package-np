## This is the serial version of npreglcaic_npRmpi.R for comparison
## purposes (bandwidth ought to be identical, timing may
## differ). Study the differences between this file and its MPI
## counterpart for insight about your own problems.

library(np)
options(np.messages=FALSE)

set.seed(42)

n <- 5000

x <- runif(n)
z1 <- rbinom(n,1,.5)
z2 <- rbinom(n,1,.5)
y <- cos(2*pi*x) + z1 + rnorm(n,sd=.25)
z1 <- factor(z1)
z2 <- factor(z2)
mydat <- data.frame(y,x,z1,z2)
rm(x,y,z1,z2)

## A regression example (local constant, aic cross-validation)  

t <- system.time(bw <- npregbw(y~x+z1+z2,
                               regtype="lc",
                               bwmethod="cv.aic",
                               data=mydat))

summary(bw)

t <- t + system.time(model <- npreg(bws=bw,
                                    data=mydat))

summary(model)

cat("Elapsed time =", t[3], "\n")

