## This is the serial version of npplreg_npRmpi.R for comparison
## purposes (bandwidth ought to be identical, timing may
## differ). Study the differences between this file and its MPI
## counterpart for insight about your own problems.

library(np)
options(np.messages=FALSE)

## Generate some data

set.seed(42)

n <- 1000

x1 <- rnorm(n)
x2 <- rbinom(n, 5, .3)

z1 <- rbinom(n, 2, .3)
z2 <- rnorm(n)

y <- 1 + x1 + x2 + z1 + sin(z2) + rnorm(n)

x2 <- ordered(x2)
z1 <- ordered(z1)

## Partially linear model

t <- system.time(bw <- npplregbw(formula=y~x1+x2|z1+z2))

summary(bw)

t <- t + system.time(pl <- npplreg(bws=bw))

summary(pl)

cat("Elapsed time =", t[3], "\n")

