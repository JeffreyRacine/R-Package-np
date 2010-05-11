sigtest <- function(In,
                    In.bootstrap,
                    P,
                    bws = NA,
                    ixvar,
                    boot.method,
                    boot.type,
                    boot.num){

  tsig <- list(In = In,
               In.bootstrap = In.bootstrap,
               P = P,
               bws = bws,
               ixvar = ixvar,
               pmethod = switch(boot.method,
                 "iid" = "IID",
                 "wild" ="Wild",
                 "wild-rademacher" = "Rademacher Wild"),
               ptype = boot.type,
               boot.num = boot.num)
  
  tsig$reject <- rep('', length(In))
  tsig$rejectNum <- rep(NA, length(In))

  tsig$reject[a <- (P < 0.1)] <- '.'
  tsig$rejectNum[a] <- 10

  tsig$reject[a <- (P < 0.05)] <- '*'
  tsig$rejectNum[a] <- 5

  tsig$reject[a <- (P < 0.01)] <- '**'
  tsig$rejectNum[a] <- 1

  tsig$reject[a <- (P < 0.001)] <- '***'
  tsig$rejectNum[a] <- 0.1

  class(tsig) = "sigtest"

  tsig
}

print.sigtest <- function(x, ...){
  cat("\nKernel Regression Significance Test",
      "\nType ", x$ptype," Test with ", x$pmethod, " Bootstrap (",x$boot.num,
      " replications)",
      "\nExplanatory variables tested for significance:\n",
      paste(paste(x$bws$xnames[x$ixvar]," (",x$ixvar,")", sep=""), collapse=", "),"\n\n",
      sep="")
      

  if (!identical(x$bws,NA))
    print(matrix(x$bws$bw,ncol=x$bws$ndim,dimnames=list(paste(x$bws$pscaling,":",sep=""),
                                            x$bws$xnames)))

  maxNameLen <- max(nc <- nchar(nm <- x$bws$xnames[x$ixvar]))

  cat("\nSignificance Tests\n")
  cat("P Value:", paste("\n", nm, ' ', blank(maxNameLen-nc), format.pval(x$P),
                        " ", x$reject, sep=''))
  cat("\n---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
}

summary.sigtest <- function(object, ...) {
  print(object)
}
