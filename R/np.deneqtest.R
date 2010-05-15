## Function that implements the multivariate density equality test
## described in Li, Q., E. Maasoumi, and J.S. Racine (2009), “A
## Nonparametric Test for Equality of Distributions with Mixed
## Categorical and Continuous Data,” Journal of Econometrics, Volume
## 148, pp 186-200.

npdeneqtest <- function(x = NULL,
                        y = NULL,
                        bw.x = NULL,
                        bw.y = NULL,
                        boot.num = 399,
                        random.seed = 42,
                        ...) {

  ## Some testing of input values

  if(is.null(x) || is.null(y)) stop(" you must provide x and y data")
  if(!is.data.frame(x) || !is.data.frame(y)) stop(" x and y must be data frames")
  if(!identical(names(data.frame(x)),names(data.frame(y)))) stop(" data frames x and y must have identical variable names")
  if(boot.num < 9) stop(" number of bootstrap replications must be >= 9")

  if(is.null(bw.x) || is.null(bw.y)) {
    bw.x <- npudensbw(dat=x,...)
    bw.y <- npudensbw(dat=y,...)     
  }

  ## Save seed prior to setting

  if(exists(".Random.seed", .GlobalEnv)) {
    save.seed <- get(".Random.seed", .GlobalEnv)
    exists.seed = TRUE
  } else {
    exists.seed = FALSE
  }

  set.seed(random.seed)

  ## First, define test statistic function. This will return the
  ## standardized and unstandardized test statistic along with its
  ## estimated variance.

  teststat <- function(x,y,bw.x,bw.y) {

    ## Get n1 and n2, number of rows in x and y

    n1 <- nrow(x)
    n2 <- nrow(y)

    ## First, compute the In statistic

    sum.1 <- sum(npksum(txdat=x,
                        bws=bw.x,
                        leave.one.out=TRUE,
                        bandwidth.divide=TRUE)$ksum)
    
    sum.2 <- sum(npksum(txdat=y,
                        bws=bw.y,
                        leave.one.out=TRUE,
                        bandwidth.divide=TRUE)$ksum)
    
    sum.3 <- sum(npksum(txdat=x,
                        exdat=y,
                        bws=bw.x,
                        leave.one.out=FALSE,
                        bandwidth.divide=TRUE)$ksum)

    ## sum.4 and sum.3 are identical...
    
    In <- sum.1/(n1*(n1-1))+sum.2/(n2*(n2-1))-2*sum.3/(n1*n2)

    ## Next, compute sigma^2_n

    sum.1 <- sum(npksum(txdat=x,
                        bws=bw.x,
                        kernel.pow=2,
                        leave.one.out=TRUE,
                        bandwidth.divide=TRUE)$ksum)
    
    sum.2 <- sum(npksum(txdat=y,
                        bws=bw.y,
                        kernel.pow=2,
                        leave.one.out=TRUE,
                        bandwidth.divide=TRUE)$ksum)
    
    sum.3 <- sum(npksum(txdat=x,
                        exdat=y,
                        bws=bw.x,
                        kernel.pow=2,
                        leave.one.out=FALSE,
                        bandwidth.divide=TRUE)$ksum)
    
    ## sum.4 and sum.3 are identical

    sigma2.n<- 2*(sum.1/(n1^2*(n1-1)^2)+sum.2/(n2^2*(n2-1)^2)+2*sum.3/(n1^2*n2^2))

    ## Finally, compute Tn, the standardized statistic

    Tn <- In/sqrt(sigma2.n)
    
    return(list(Tn=Tn,In=In))
    
  } ## End of test statistic

  ## Now write a bootstrap function for the test statistic
  
  teststat.boot <- function(x,y,bw.x,bw.y) {
    n1 <- nrow(x)
    n2 <- nrow(y)
    ## Resample from pooled data
    z <- data.frame(rbind(x,y))
    x.bootstrap <- data.frame(z[sample(nrow(z),size=n1,replace=TRUE),])      
    y.bootstrap <- data.frame(z[sample(nrow(z),size=n2,replace=TRUE),])      
    output.boot <- teststat(x.bootstrap,y.bootstrap,bw.x,bw.y)
    return(list(Tn=output.boot$Tn,
                In=output.boot$In))
  }
  
  Tn.vector <- numeric(boot.num)
  In.vector <- numeric(boot.num)
  
  console <- newLineConsole()

  for(i in 1:boot.num) {
    console <- printClear(console)
    console <- printPush(paste(sep="", "Bootstrap replication ",
                               i, "/", boot.num, "..."), console)
    output.boot <- teststat.boot(x,y,bw.x,bw.y)
    Tn.vector[i] <- output.boot$Tn
    In.vector[i] <- output.boot$In
  }
  
  console <- printClear(console)
  console <- printPop(console)  

  ## Compute the test statistic
  
  output <- teststat(x,y,bw.x,bw.y)
  
  ## Compute empirical P-values - the number of resampled statistics
  ## more extreme than the original statistic
  
  Tn.P <- mean(ifelse(Tn.vector>output$Tn,1,0))
  In.P <- mean(ifelse(In.vector>output$In,1,0))
  
  ## Restore seed

  if(exists.seed) assign(".Random.seed", save.seed, .GlobalEnv)
  
  deneqtest(Tn=output$Tn,
            In=output$In,
            Tn.bootstrap=Tn.vector,
            In.bootstrap=In.vector,                
            Tn.P=Tn.P,
            In.P=In.P,
            boot.num=boot.num)
  
}
