## This is the serial version of npcmstest_npRmpi.R for comparison
## purposes (bandwidth ought to be identical, timing may
## differ). Study the differences between this file and its MPI
## counterpart for insight about your own problems.

library(np)
options(np.messages=FALSE)

data("oecdpanel")
attach(oecdpanel)
n <- nrow(oecdpanel)

oecd <- factor(oecd)
year <- factor(year)

model <- lm(growth ~ oecd +
            year +
            initgdp +
            I(initgdp^2) +
            I(initgdp^3) +
            I(initgdp^4) +
            popgro +
            inv +
            humancap +
            I(humancap^2) +
            I(humancap^3) - 1, 
            x=TRUE, 
            y=TRUE)

X <- data.frame(oecd, year, initgdp, popgro, inv, humancap)

## Consistent model specification test (we override defaults for
## demonstration purposes - don't do this for real problems).

t <- system.time(output <- npcmstest(model = model,
                                     xdat = X,
                                     ydat = growth,
                                     nmulti=1,
                                     ftol=.01,
                                     tol=.01))

output

cat("Elapsed time =", t[3], "\n")
