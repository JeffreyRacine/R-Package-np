## Function that implements the multivariate density equality test
## described in Li, Q., E. Maasoumi, and J.S. Racine (2009), "A
## Nonparametric Test for Equality of Distributions with Mixed
## Categorical and Continuous Data," Journal of Econometrics, Volume
## 148, pp 186-200.

.npdeneq_count_plan <- function(index, pool.n) {
  counts <- tabulate(index, nbins = pool.n)
  support <- which(counts > 0L)
  list(index = support, counts = as.numeric(counts[support]))
}

.npdeneq_count_compression_eligible <- function(bw) {
  if (is.numeric(bw) || is.integer(bw))
    return(TRUE)

  kbw <- tryCatch(kbandwidth(bw), error = function(e) NULL)
  if (is.null(kbw) ||
      !identical(kbw[["type", exact = TRUE]], "fixed"))
    return(FALSE)

  lower <- kbw[["ckerlb", exact = TRUE]]
  upper <- kbw[["ckerub", exact = TRUE]]
  !any(is.finite(lower) | is.finite(upper))
}

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

  .np_progress_note("Computing bandwidths")

  if(is.null(bw.x) || is.null(bw.y)) {
    bw.x <- .np_progress_with_legacy_suppressed(npudensbw(dat=x,...))
    bw.y <- .np_progress_with_legacy_suppressed(npudensbw(dat=y,...))
  }

  ## Save seed prior to setting

  seed.state <- .np_seed_enter(random.seed)
  on.exit(.np_seed_exit(seed.state, remove_if_absent = TRUE), add = TRUE)


  ## First, define test statistic function. This will return the
  ## standardized and unstandardized test statistic along with its
  ## estimated variance.

  teststat <- function(x,y,bw.x,bw.y) {

    ## Get n1 and n2, number of rows in x and y

    n1 <- nrow(x)
    n2 <- nrow(y)

    ## First, compute the In statistic

    ksum.1 <- .npksum_power12(txdat=x,
                              bws=bw.x,
                              leave.one.out=TRUE,
                              bandwidth.divide=TRUE)
    sum.1 <- sum(ksum.1$ksum)
    sum2.1 <- sum(ksum.1$ksum.power2)

    ksum.2 <- .npksum_power12(txdat=y,
                              bws=bw.y,
                              leave.one.out=TRUE,
                              bandwidth.divide=TRUE)
    sum.2 <- sum(ksum.2$ksum)
    sum2.2 <- sum(ksum.2$ksum.power2)

    ksum.3 <- .npksum_power12(txdat=x,
                              exdat=y,
                              bws=bw.x,
                              leave.one.out=FALSE,
                              bandwidth.divide=TRUE)
    sum.3 <- sum(ksum.3$ksum)
    sum2.3 <- sum(ksum.3$ksum.power2)

    ## sum.4 and sum.3 are identical...
    
    In <- sum.1/(n1*(n1-1))+sum.2/(n2*(n2-1))-2*sum.3/(n1*n2)

    ## Next, compute sigma^2_n

    ## sum.4 and sum.3 are identical

    sigma2.n<- 2*(sum2.1/(n1^2*(n1-1)^2)+sum2.2/(n2^2*(n2-1)^2)+2*sum2.3/(n1^2*n2^2))

    ## Finally, compute Tn, the standardized statistic

    Tn <- In/sqrt(sigma2.n)
    
    return(list(Tn=Tn,In=In))
    
  } ## End of test statistic

  teststat.counted <- function(z, x.count, y.count, bw.x, bw.y, n1, n2) {
    x.support <- z[x.count$index, , drop = FALSE]
    y.support <- z[y.count$index, , drop = FALSE]

    ksum.1 <- .npksum_power12_weighted(
      txdat = x.support,
      counts = x.count$counts,
      bws = bw.x,
      bandwidth.divide = TRUE
    )
    sum.1 <- sum(x.count$counts *
                 (as.numeric(ksum.1$ksum) -
                  as.numeric(self.diagonal.x$ksum)))
    sum2.1 <- sum(x.count$counts *
                  (as.numeric(ksum.1$ksum.power2) -
                   as.numeric(self.diagonal.x$ksum.power2)))

    ksum.2 <- .npksum_power12_weighted(
      txdat = y.support,
      counts = y.count$counts,
      bws = bw.y,
      bandwidth.divide = TRUE
    )
    sum.2 <- sum(y.count$counts *
                 (as.numeric(ksum.2$ksum) -
                  as.numeric(self.diagonal.y$ksum)))
    sum2.2 <- sum(y.count$counts *
                  (as.numeric(ksum.2$ksum.power2) -
                   as.numeric(self.diagonal.y$ksum.power2)))

    ksum.3 <- .npksum_power12_weighted(
      txdat = x.support,
      counts = x.count$counts,
      exdat = y.support,
      bws = bw.x,
      bandwidth.divide = TRUE
    )
    sum.3 <- sum(y.count$counts * as.numeric(ksum.3$ksum))
    sum2.3 <- sum(y.count$counts * as.numeric(ksum.3$ksum.power2))

    In <- sum.1 / (n1 * (n1 - 1)) +
      sum.2 / (n2 * (n2 - 1)) -
      2 * sum.3 / (n1 * n2)
    sigma2.n <- 2 * (
      sum2.1 / (n1^2 * (n1 - 1)^2) +
      sum2.2 / (n2^2 * (n2 - 1)^2) +
      2 * sum2.3 / (n1^2 * n2^2)
    )

    list(Tn = In / sqrt(sigma2.n), In = In)
  }

  compress.bootstrap <-
    .npdeneq_count_compression_eligible(bw.x) &&
    .npdeneq_count_compression_eligible(bw.y)
  bootstrap.pool <- data.frame(rbind(x, y))
  self.diagonal.x <- self.diagonal.y <- NULL
  if (compress.bootstrap) {
    self.diagonal.x <- .npksum_power12(
      txdat = bootstrap.pool[1L, , drop = FALSE],
      bws = bw.x,
      bandwidth.divide = TRUE
    )
    self.diagonal.y <- .npksum_power12(
      txdat = bootstrap.pool[1L, , drop = FALSE],
      bws = bw.y,
      bandwidth.divide = TRUE
    )
  }

  ## Now write a bootstrap function for the test statistic
  
  teststat.boot <- function(x,y,bw.x,bw.y) {
    n1 <- nrow(x)
    n2 <- nrow(y)
    ## Resample from pooled data
    z <- bootstrap.pool
    x.index <- sample.int(nrow(z), size = n1, replace = TRUE)
    y.index <- sample.int(nrow(z), size = n2, replace = TRUE)
    output.boot <- if (compress.bootstrap) {
      teststat.counted(
        z = z,
        x.count = .npdeneq_count_plan(x.index, nrow(z)),
        y.count = .npdeneq_count_plan(y.index, nrow(z)),
        bw.x = bw.x,
        bw.y = bw.y,
        n1 = n1,
        n2 = n2
      )
    } else {
      x.bootstrap <- data.frame(z[x.index, , drop = FALSE])
      y.bootstrap <- data.frame(z[y.index, , drop = FALSE])
      teststat(x.bootstrap, y.bootstrap, bw.x, bw.y)
    }
    return(list(Tn=output.boot$Tn,
                In=output.boot$In))
  }
  
  Tn.vector <- numeric(boot.num)
  In.vector <- numeric(boot.num)

  progress <- .np_progress_begin("Bootstrap replications", total = boot.num, surface = "bootstrap")

  for (i in seq_len(boot.num)) {
    output.boot <- teststat.boot(x,y,bw.x,bw.y)
    Tn.vector[i] <- output.boot$Tn
    In.vector[i] <- output.boot$In
    progress <- .np_progress_step(progress, done = i)
  }

  progress <- .np_progress_end(progress)

  ## Compute the test statistic
  
  output <- teststat(x,y,bw.x,bw.y)
  
  ## Compute empirical P-values - the number of resampled statistics
  ## more extreme than the original statistic
  
  Tn.P <- mean(Tn.vector > output$Tn)
  In.P <- mean(In.vector > output$In)
  
  ## Restore seed

  .np_seed_exit(seed.state, remove_if_absent = TRUE)
  
  deneqtest(Tn=output$Tn,
            In=output$In,
            Tn.bootstrap=Tn.vector,
            In.bootstrap=In.vector,                
            Tn.P=Tn.P,
            In.P=In.P,
            boot.num=boot.num)
  
}
