## This is the serial version of npindexks_npRmpi.R for comparison
## purposes (bandwidth ought to be identical, timing may
## differ). Study the differences between this file and its MPI
## counterpart for insight about your own problems.

library(np)
options(np.messages=FALSE)

## Generate some data

set.seed(42)

n <- 5000

x <- rchisq(n, df=3)
x1 <- (ifelse(x < 6, x, 6) - 2.348)/1.511
x <- rnorm(n)
x2 <- ifelse(abs(x) < 2 , x, 2) / 0.8796
y <- ifelse(x1 + x2 + rnorm(n) > 0, 1, 0)
mydat <- data.frame(x1,x2,y)
rm(x,x1,x2,y)
     
## A single index model example (Klein & Spady, binary y)

t <- system.time(bw <- npindexbw(formula=y~x1+x2,
                                 method="kleinspady",
                                 data=mydat))

summary(bw)

t <- t + system.time(model <- npindex(bws=bw, gradients=TRUE))

summary(model)

cat("Elapsed time =", t[3], "\n")

