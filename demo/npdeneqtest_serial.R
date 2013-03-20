## This is the serial version of npdeneqtest_npRmpi.R for comparison
## purposes (bandwidth ought to be identical, timing may
## differ). Study the differences between this file and its MPI
## counterpart for insight about your own problems.

library(np)
options(np.messages=FALSE)

## Generate some data

set.seed(42)

n <- 2500

sample.A <- data.frame(x=rnorm(n))
sample.B <- data.frame(x=rnorm(n))

## A consistent density equality test example

t <- system.time(output <- npdeneqtest(sample.A,sample.B,boot.num=99))

output

cat("Elapsed time =", t[3], "\n")
