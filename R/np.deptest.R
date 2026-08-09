## Function that implements the bivariate dependence metric described
## in Maasoumi, E.~and J.S.~Racine (2002), "Entropy and Predictability
## of Stock Market Returns," Journal of Econometrics, March, Volume
## 107, Issue 2, pp 291-312.

npdeptest <- function(data.x = NULL,
                      data.y = NULL,
                      method=c("integration","summation"),
                      bootstrap = TRUE,
                      boot.num = 399,
                      random.seed = 42) {
  
  ## Trap fatal errors

  if(is.data.frame(data.x)||is.data.frame(data.y)) stop(" you must enter two data vectors (and not data frames)")
  if(is.factor(data.x)||is.factor(data.y)) stop(" does not support factors")
  if(is.null(data.x)||is.null(data.y)) stop(" you must enter x and y data vectors")
  if(ncol(data.frame(data.x)) != 1) stop(" data must have one dimension only")
  if(length(data.x)!=length(data.y)) stop(" data vectors must be of equal length")
  if(boot.num < 9) stop(" number of bootstrap replications must be >= 9")

  method <- match.arg(method)

  ## Save seed prior to setting

  seed.state <- .np_seed_enter(random.seed)
  on.exit(.np_seed_exit(seed.state, remove_if_absent = TRUE), add = TRUE)

  .np_progress_note("Computing bandwidths")

  ## If the variable is a time series convert to type numeric

  if(is.ts(data.x)) data.x <- as.numeric(data.x)
  if(is.ts(data.y)) data.y <- as.numeric(data.y)

  ## Remove any NAs from paired data

  tmp <- na.omit(data.frame(data.x,data.y))
  data.x <- tmp$data.x
  data.y <- tmp$data.y
  rm(tmp)

  ## Define the metric entropy function

  Srho.bivar <- function(x.dat = NULL,
                         y.dat = NULL,
                         bw.x = NULL,
                         bw.y = NULL,
                         bw.joint = NULL,
                         method=c("integration","summation")) {
    
    if(is.null(bw.x)||is.null(bw.y)||is.null(bw.joint)) stop(" you must provide numeric bandwidths for f(x), f(y) and f(x,y)")
    
    method <- match.arg(method)
    
    if(method=="summation") {
      
      ## Summation version: \sum_i(1-sqrt(f(x_i)f(y_i)/f(x_i,y_i)))^2.
      ## The fixed Gaussian marginal and joint sums share one native,
      ## constant-workspace traversal.
      
      return(.np_entropy_bivariate_gaussian_summation(
        x.dat, y.dat, bw.x, bw.y, bw.joint
      ))
      
    } else {
      
      ## Integration version:
      ## \int\int (sqrt(f(x,y))-sqrt(f(x))sqrt(f(y)))^2 dx dy
      
      domain <- .np_entropy_bivariate_domain(
        x.dat, y.dat, bw.x, bw.y, bw.joint
      )
      
      return(.np_entropy_bivariate_integral(
        x.dat = x.dat,
        y.dat = y.dat,
        bw.x = bw.x,
        bw.y = bw.y,
        bw.joint = bw.joint,
        lower = domain$lower,
        upper = domain$upper
      ))
      
    }

  } ## end of Srho.bivar function
  
  ## Compute and save bandwidths (save for bootstrapping if requested)

  bw.data.x <- .np_progress_with_legacy_suppressed(npudensbw(~data.x))$bw
  bw.data.y <- .np_progress_with_legacy_suppressed(npudensbw(~data.y))$bw
  bw.joint <- .np_progress_with_legacy_suppressed(npudensbw(~data.x+data.y))$bw

  .np_progress_note("Constructing metric entropy")
  
  Srho.vec <- Srho.bivar(data.x,data.y,bw.data.x,bw.data.y,bw.joint,method=method)

  ## Bootstrap if requested - null is independence so simple iid
  ## index resampling under replacement is sufficient

  if(bootstrap) {

    Srho.vec.boot <- numeric(boot.num)
    progress <- .np_progress_begin("Bootstrap replications", total = boot.num, surface = "bootstrap")

    if (method == "summation") {
      sample.size <- length(data.x)
      chunk.size <- .np_entropy_count_chunk_size(
        sample.size, bytes.per.support = 4, max.chunk = 16L
      )
      for (start in seq.int(1L, boot.num, by = chunk.size)) {
        index <- start:min(boot.num, start + chunk.size - 1L)
        bootstrap.index <- matrix(
          NA_integer_, nrow = length(index), ncol = sample.size
        )
        for (row in seq_along(index)) {
          bootstrap.index[row, ] <- sample.int(
            sample.size, replace = TRUE
          )
        }
        Srho.vec.boot[index] <-
          .np_entropy_bivariate_gaussian_summation_xindex(
            data.x, data.y, bootstrap.index,
            bw.data.x, bw.data.y, bw.joint
          )
        for (done in index)
          progress <- .np_progress_step(progress, done = done)
      }
    } else {
      for (b in seq_len(boot.num)) {
        ## Break systematic relationship between x and y (null)

        data.x.boot <- data.x[sample.int(length(data.x), replace = TRUE)]

        Srho.vec.boot[b] <- Srho.bivar(data.x.boot,data.y,bw.data.x,bw.data.y,bw.joint,method=method)
        progress <- .np_progress_step(progress, done = b)
      }
    }

    progress <- .np_progress_end(progress)

    ## Compute P-values

    P <- mean(Srho.vec.boot > Srho.vec)

  }

  ## Restore seed

  .np_seed_exit(seed.state, remove_if_absent = TRUE)

  if(bootstrap) {

    deptest(Srho = Srho.vec,
            Srho.bootstrap.vec = Srho.vec.boot,
            P = P,
            bootstrap = bootstrap,
            boot.num = boot.num,
            bw.data.x = bw.data.x,
            bw.data.y = bw.data.y,
            bw.joint = bw.joint)

  } else {

    deptest(Srho = Srho.vec,
            Srho.bootstrap.vec = NULL,
            P = NULL,
            bootstrap = bootstrap,
            boot.num = NULL,
            bw.data.x = bw.data.x,
            bw.data.y = bw.data.y,
            bw.joint = bw.joint)

  }

}
