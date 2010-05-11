deneqtest <- function(Tn,
                      In,
                      Tn.bootstrap,
                      In.bootstrap,
                      Tn.P,
                      In.P,                      
                      boot.num) {

  tdeneq = list(Tn=Tn,
    In=In,
    Tn.bootstrap=Tn.bootstrap,
    In.bootstrap=In.bootstrap,                
    Tn.P=Tn.P,
    In.P=In.P,    
    boot.num=boot.num)

  reject <- ''
  
  if (Tn.P < 0.1)
    reject <- '.'

  if (Tn.P < 0.05)
    reject <- '*'

  if (Tn.P < 0.01)
    reject <- '**'

  if (Tn.P < 0.001)
    reject <- '***'

  tdeneq$reject <- reject
  tdeneq$rejectNum <- switch(reject,
                           '.' = 10,
                           '*' = 5,
                           '**' = 1,
                           '***' = 0.1)
  
  class(tdeneq) <- "deneqtest"

  tdeneq
}

print.deneqtest <- function(x, ...){
  cat("\nConsistent Density Equality Test",
      paste("\n", format(x$boot.num), " Bootstrap Replications",sep=""),
      "\n\nTest Statistic ", sQuote("Tn"), ": ",
      format(x$Tn), "\tP Value: ", format.pval(x$Tn.P)," ", x$reject,
      "\n---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1",
      ifelse(x$reject == '', "\nFail to reject the null of equality at the 10% level",
             paste("\nNull of equality is rejected at the ", x$rejectNum, "% level", sep="")),
      "\n\n", sep="")
}

summary.deneqtest <- function(object, ...){
  print(object)
}
