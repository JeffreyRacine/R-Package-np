cmstest <- function(Jn, In, Omega.hat, sd,
                    q.90, q.95, q.99, P,
                    bws,
                    distribution,
                    Jn.bootstrap = NA, In.bootstrap = NA,
                    pivot,
                    model,
                    boot.method,
                    boot.num,
                    na.index = NULL){

  tcms = list(Jn = Jn,
    In = In,
    Omega.hat = Omega.hat,
    sd = sd,
    q.90 = q.90,
    q.95 = q.95,
    q.99 = q.99,
    P = P,
    bws = bws,
    Jn.bootstrap = Jn.bootstrap,
    In.bootstrap = In.bootstrap,
    pivot = pivot,
    pdistribution = switch(distribution,
      bootstrap = "Bootstrap",
      asymptotic = "Asymptotic"),
    pcall = paste(deparse(model$call),collapse=""),
    pmethod = switch(boot.method,
      "iid" = "IID",
      "wild" ="Wild",
      "wild-rademacher" = "Rademacher Wild"),
    boot.num = boot.num,
    na.index = na.index)

  ##Sn = if(pivot) Jn else In

  reject <- ''
  
  if (P < 0.1)
    reject <- '.'

  if (P < 0.05)
    reject <- '*'

  if (P < 0.01)
    reject <- '**'

  if (P < 0.001)
    reject <- '***'

  tcms$reject <- reject
  tcms$rejectNum <- switch(reject,
                           '.' = 10,
                           '*' = 5,
                           '**' = 1,
                           '***' = 0.1)
  
  class(tcms) <- "cmstest"

  tcms
}

print.cmstest <- function(x, ...){
  cat("\nConsistent Model Specification Test",
      "\nParametric null model: ", x$pcall,
      "\nNumber of regressors: ", length(x$bws$bw),"\n",
      if(x$pdistribution == "Bootstrap"){
        paste(x$pmethod, " Bootstrap ",
              "(", x$boot.num, " replications)",sep="")
      } else paste(x$pdistribution, "Distribution"),
      "\n\nTest Statistic ",sQuote(ifelse(x$pivot,"Jn","In")), ": ",
      format(ifelse(x$pivot, x$Jn, x$In)), "\tP Value: ", format.pval(x$P)," ", x$reject,
      "\n---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1",
      ifelse(x$reject == '', "\nFail to reject the null of correct specification at the 10% level",
             paste("\nNull of correct specification is rejected at the ", x$rejectNum, "% level", sep="")),
      "\n\n", sep="")
}

summary.cmstest <- function(object, ...){
  print(object)
}
