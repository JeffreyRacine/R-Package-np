## This is the serial version of npconmode_npRmpi.R for comparison
## purposes (bandwidth ought to be identical, timing may
## differ). Study the differences between this file and its MPI
## counterpart for insight about your own problems.

library(np)
options(np.messages=FALSE)

library(MASS)
data(birthwt)
n <- nrow(birthwt)

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
