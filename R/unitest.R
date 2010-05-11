unitest <- function(Srho,
                       Srho.bootstrap,
                       P,
                       bootstrap,
                       boot.num,
                       bw.x,
                       bw.y) {
  
  tuni = list(Srho = Srho,
    Srho.bootstrap = Srho.bootstrap,
    P = P,
    bootstrap = bootstrap,
    boot.num = boot.num,
    bw.x = bw.x,
    bw.y = bw.y)

  if(bootstrap) {
  
    reject <- ''
    
    if (P < 0.1)
      reject <- '.'
    
    if (P < 0.05)
      reject <- '*'
    
    if (P < 0.01)
      reject <- '**'
    
    if (P < 0.001)
      reject <- '***'
    
    tuni$reject <- reject
    tuni$rejectNum <- switch(reject,
                             '.' = 10,
                             '*' = 5,
                             '**' = 1,
                             '***' = 0.1)

  }
  
  class(tuni) <- "unitest"

  tuni

}

print.unitest <- function(x, ...){
  if(x$bootstrap) {
    cat("\nConsistent Univariate Entropy Density Equality Test",
        paste("\n", format(x$boot.num), " Bootstrap Replications",sep=""),
        "\n\nTest Statistic ", sQuote("Srho"), ": ",
        format(x$Srho), "\tP Value: ", format.pval(x$P)," ", x$reject,
        "\n---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1",
        ifelse(x$reject == '', "\nFail to reject the null of equality at the 10% level",
               paste("\nNull of equality is rejected at the ", x$rejectNum, "% level", sep="")),
        "\n\n", sep="")
  } else {
    cat("\nConsistent Univariate Density Difference Metric Entropy",
        "\n\nMetric Entropy ", sQuote("Srho"),": ",
        format(x$Srho),"\n\n", sep="")
  }
}

summary.unitest <- function(object, ...){
  print(object)
}
