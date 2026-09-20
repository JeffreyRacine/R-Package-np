# Test-only arithmetic: independently reconstruct each admitted NN sample.
# Columns are held-out occurrences; rows are donors. No package kernel/NN
# builder or leave-one-out implementation participates in this oracle.
workhorse_nn_loo_weights <- function(training, bandwidth, bwtype,
                                      kernel="beta", order=2L) {
  n <- nrow(training)
  result <- matrix(0,n,n)
  for(i in seq_len(n)) {
    donors <- setdiff(seq_len(n),i)
    weight <- rep(1,length(donors))
    for(d in seq_len(ncol(training))) {
      x <- training[[d]]
      h <- if(bwtype=="fixed") rep(bandwidth[d],length(donors)) else
        if(bwtype=="generalized_nn")
          rep(sort(abs(x[donors]-x[i]))[bandwidth[d]],length(donors)) else
          vapply(donors,function(j)
            sort(abs(x[setdiff(donors,j)]-x[j]))[bandwidth[d]],numeric(1))
      if(kernel=="beta") {
        m <- order/2L
        value <- rep(0,length(donors))
        for(s in seq_len(m))
          value <- value+(-1)^(s+1)*choose(m,s)*
            dbeta(x[donors],1+x[i]/(s*h^2),1+(1-x[i])/(s*h^2))
      } else {
        stopifnot(kernel=="gaussian")
        u <- (x[donors]-x[i])/h
        polynomial <- switch(as.character(order),`2`=1,`4`=1.5-.5*u^2,
          `6`=1.875-1.25*u^2+.125*u^4,
          # Retain the incumbent eighth-order polynomial's decimal
          # coefficient; this test changes NN support, not kernel constants.
          `8`=35/16-35/16*u^2+7/16*u^4-0.02083333333*u^6)
        value <- dnorm(u)*polynomial/h
      }
      weight <- weight*value
    }
    result[donors,i] <- weight
  }
  result
}
