symtest <- function(Srho, Srho.bootstrap, P, boot.num, data.rotate, bw) {
  
  tsym = list(Srho = Srho,
    Srho.bootstrap = Srho.bootstrap,
    P = P,
    boot.num = boot.num,
    data.rotate = data.rotate,
    bw = bw)
  
  reject <- ''
  
  if (P < 0.1)
    reject <- '.'

  if (P < 0.05)
    reject <- '*'

  if (P < 0.01)
    reject <- '**'

  if (P < 0.001)
    reject <- '***'

  tsym$reject <- reject
  tsym$rejectNum <- switch(reject,
                           '.' = 10,
                           '*' = 5,
                           '**' = 1,
                           '***' = 0.1)
  
  class(tsym) <- "symtest"

  tsym
}

print.symtest <- function(x, ...){
  cat("\nConsistent Entropy Asymmetry Test",
      paste("\n", format(x$boot.num), " Bootstrap Replications",sep=""),
      "\n\nTest Statistic ", sQuote("Srho"), ": ",
      format(x$Srho), "\tP Value: ", format.pval(x$P)," ", x$reject,
      "\n---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1",
      ifelse(x$reject == '', "\nFail to reject the null of symmetry at the 10% level",
             paste("\nNull of symmetry is rejected at the ", x$rejectNum, "% level", sep="")),
      "\n\n", sep="")
}

summary.symtest <- function(object, ...){
  print(object)
}
