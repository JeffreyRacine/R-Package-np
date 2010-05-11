deptest <- function(Srho,
                    Srho.bootstrap.vec,
                    P,
                    bootstrap,
                    boot.num,
                    bw.data.x,
                    bw.data.y,
                    bw.joint) {
  
  tdep = list(Srho = Srho,
    Srho.bootstrap.vec = Srho.bootstrap.vec,
    P = P,
    bootstrap = bootstrap,
    boot.num = boot.num,
    bw.data.x = bw.data.x,
    bw.data.y = bw.data.y,
    bw.joint = bw.joint)

  if(is.null(P)) {

    tdep$reject <- NA
    tdep$rejectNum <- NA

  } else {

    reject <- ''
    
    if (P < 0.1)
      reject <- '.'
    
    if (P < 0.05)
      reject <- '*'
    
    if (P < 0.01)
      reject <- '**'
    
    if (P < 0.001)
      reject <- '***'
    
    tdep$reject <- reject
    tdep$rejectNum <- switch(reject,
                             '.' = 10,
                             '*' = 5,
                             '**' = 1,
                             '***' = 0.1)
    
  }
  
  class(tdep) <- "deptest"

  tdep
}

print.deptest <- function(x, ...){
  if(x$bootstrap) {
    cat("\nConsistent Metric Entropy Test for Dependence",
        paste("\n", format(x$boot.num), " Bootstrap Replications",sep=""))
    cat(paste("\n\nTest Statistic ", sQuote("Srho"), ": ",
              format(x$Srho), "\tP Value: ", format.pval(x$P)," ", x$reject,sep=""))
    cat("\n---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
    cat(ifelse(x$reject == '', paste("\nFail to reject the null of independence at the 10% level\n\n",sep=""),
               paste("\nNull of independence is rejected at the ", x$rejectNum, "% level\n\n", sep="")))
  } else {
    cat("\nConsistent Dependence Metric Entropy",
        paste("\nStatistic ", sQuote("Srho"),sep=""), ": ",
        format(x$Srho),"\n\n",sep="")
  }
}

summary.deptest <- function(object, ...){
  print(object)
}
