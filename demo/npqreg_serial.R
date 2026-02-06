## This is the serial version of npqreg_npRmpi.R for comparison
## purposes (bandwidth ought to be identical, timing may
## differ). Study the differences between this file and its MPI
## counterpart for insight about your own problems.

library(np)
options(np.messages=FALSE)

data("Italy")

n <- 1008 # n <- nrow(Italy)

## A quantile regression example

t <- system.time(bw <- npcdistbw(gdp~ordered(year),data=Italy))

summary(bw)

t <- t + system.time(model.q0.25 <- npqreg(bws=bw, tau=0.25))
t <- t + system.time(model.q0.50 <- npqreg(bws=bw, tau=0.50))
t <- t + system.time(model.q0.75 <- npqreg(bws=bw, tau=0.75))

cat("Elapsed time =", t[3], "\n")
