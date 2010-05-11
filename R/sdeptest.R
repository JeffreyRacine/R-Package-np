sdeptest <- function(Srho,
                    Srho.cumulant,
                    Srho.bootstrap.mat,
                    Srho.cumulant.bootstrap.mat,
                    P,
                    P.cumulant,
                    bootstrap,
                    boot.num,
                    lag.num,
                    bw.y,
                    bw.y.lag,
                    bw.joint) {
  
  tsdep = list(Srho = Srho,
    Srho.cumulant = Srho.cumulant,
    Srho.bootstrap.mat = Srho.bootstrap.mat,
    Srho.cumulant.bootstrap.mat = Srho.cumulant.bootstrap.mat,
    P = P,
    P.cumulant = P.cumulant,
    bootstrap = bootstrap,
    boot.num = boot.num,
    lag.num = lag.num,
    bw.y = bw.y,
    bw.y.lag = bw.y.lag,
    bw.joint = bw.joint)

  if(tsdep$bootstrap) {

    tsdep$reject <- numeric()
    tsdep$rejectNum <- character()

    for(k in 1:lag.num) {
      
      reject <- ''
      
      if (P[k] < 0.1)
        reject <- '.'
      
      if (P[k] < 0.05)
        reject <- '*'
      
      if (P[k] < 0.01)
        reject <- '**'
      
      if (P[k] < 0.001)
        reject <- '***'
      
      tsdep$reject[k] <- reject
      tsdep$rejectNum[k] <- switch(reject,
                                  '.' = 10,
                                  '*' = 5,
                                  '**' = 1,
                                  '***' = 0.1)

    }
    
  }
    
  class(tsdep) <- "sdeptest"

  tsdep
}

print.sdeptest <- function(x, ...){
  if(x$bootstrap) {
    cat("\nConsistent Metric Entropy Test for Nonlinear Dependence",
        paste("\n", format(x$boot.num), " Bootstrap Replications, ",
              format(x$lag.num), " Lags\n",sep=""))
    for(k in 1:x$lag.num) {
      cat(paste("\nTest Statistic ", sQuote(paste("Srho[", k , "]",sep="")), ": ",
          format(x$Srho[k]), "\tP Value: ", format.pval(x$P[k])," ", x$reject[k],sep=""))
    }

    cat("\n---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
    
    for(k in 1:x$lag.num) {
      cat(ifelse(x$reject[k] == '', paste("\nFail to reject the null of independence at lag ", k, " at the 10% level",sep=""),
                 paste("\nNull of independence is rejected at lag ", k, " at the ", x$rejectNum[k], "% level", sep="")))
    }
  } else {
    cat("\nConsistent Nonlinear Dependence Metric Entropy",
        paste("\n", format(x$lag.num), " Lags\n",sep=""))
    for(k in 1:x$lag.num) {
      cat(paste("\nStatistic ", sQuote(paste("Srho[", k , "]",sep="")), ": ",
          format(x$Srho[k]),sep=""))
    }
  }
  cat("\n\n")
}

summary.sdeptest <- function(object, ...){
  print(object)
}
