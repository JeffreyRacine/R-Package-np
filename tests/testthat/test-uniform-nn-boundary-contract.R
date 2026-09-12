local({
  # Independent dense test oracle only; production retains streamed/tree rows.
  uniform.rows <- function(train, evaluation, k, type) {
    difference <- outer(train,evaluation,"-")
    if(type == "adaptive_nn") {
      rank.distance <- abs(outer(train,train,"-"))
      diag(rank.distance) <- Inf
      radius <- apply(rank.distance,1L,function(z)sort(z)[k])
      inside <- sweep(abs(difference),1L,radius,"<")
      normalized <- sweep(difference,1L,radius,"/")
    } else {
      radius <- apply(abs(difference),2L,function(z)sort(z)[k])
      inside <- sweep(abs(difference),2L,radius,"<")
      normalized <- sweep(difference,2L,radius,"/")
    }
    list(normal=.5*inside, integral=pmin(pmax((1-normalized)/2,0),1))
  }
  run <- function() {
    old <- options(np.messages=FALSE,np.extendednn=TRUE)
    on.exit(options(old),add=TRUE)
    set.seed(42); invisible(rnorm(150)); x <- data.frame(x=runif(150))
    test_that("uniform ANN strict membership is independent of evaluation batch", {
      bw <- npudensbw(dat=x,bws=15,bwtype="adaptive_nn",ckertype="uniform",
                     bandwidth.compute=FALSE)
      distance <- abs(outer(x$x,x$x,"-")); diag(distance) <- Inf
      h <- apply(distance,1L,function(z)sort(z)[15L])
      expected <- colMeans(sweep(uniform.rows(x$x,x$x,15L,"adaptive_nn")$normal,
                                1L,h,"/"))
      for(tree in list(FALSE,TRUE,"auto")) {
        options(np.tree=tree)
        all <- fitted(npudens(bws=bw,tdat=x,edat=x))
        single <- vapply(seq_len(nrow(x)),function(j)
          fitted(npudens(bws=bw,tdat=x,edat=x[j,,drop=FALSE])),numeric(1))
        expect_equal(all,expected,tolerance=2e-13)
        expect_equal(single,expected,tolerance=2e-13)
      }
    })
    set.seed(7); x <- data.frame(x=rnorm(40))
    test_that("uniform ANN excludes the radius but keeps its immediate interior", {
      train <- data.frame(x=c(0,1,2,4,8,16,32))
      at <- data.frame(x=c(2-.Machine$double.eps,2,2+2*.Machine$double.eps))
      expected <- uniform.rows(train$x,at$x,2L,"adaptive_nn")$normal
      expect_identical(as.vector(expected[1L,]),c(.5,0,0))
      for(tree in list(FALSE,TRUE,"auto")) {
        options(np.tree=tree)
        weights <- npksum(txdat=train,exdat=at,bws=2,bwtype="adaptive_nn",
          ckertype="uniform",return.kernel.weights=TRUE)$kw
        expect_identical(unname(weights),expected)
      }
    })
    test_that("both uniform NN modes respect affine density scaling", {
      for(type in c("generalized_nn","adaptive_nn"))
        for(tree in list(FALSE,TRUE,"auto")) {
          options(np.tree=tree)
          transformed <- data.frame(x=3*x$x+1)
          bw <- npudensbw(dat=x,bws=6,bwtype=type,ckertype="uniform",
                         bandwidth.compute=FALSE)
          bz <- npudensbw(dat=transformed,bws=6,bwtype=type,ckertype="uniform",
                         bandwidth.compute=FALSE)
          expect_equal(fitted(npudens(bws=bw,tdat=x)),
                       3*fitted(npudens(bws=bz,tdat=transformed)),tolerance=2e-13)
        }
    })
    test_that("normal factors in permutation and mixed integral rows share membership", {
      i <- seq_len(48L)
      x <- data.frame(a=sin(i*.714)+i/300,b=cos(i*.417)+sin(i*.71)/3)
      at <- x[c(7L,23L,39L),]
      for(type in c("generalized_nn","adaptive_nn")) {
        a <- uniform.rows(x$a,at$a,11L,type)
        b <- uniform.rows(x$b,at$b,15L,type)
        for(tree in list(FALSE,TRUE,"auto")) {
          options(np.tree=tree)
          for(op in list(c("normal","normal"),c("normal","integral"),
                         c("integral","normal"),c("integral","integral"))) {
            z <- npksum(txdat=x,exdat=at,bws=c(11,15),bwtype=type,
              ckertype="uniform",operator=op,permutation.operator="normal",
              return.kernel.weights=TRUE,return.derivative.kernel.weights=TRUE)
            expect_equal(unname(z$kw),a[[op[1L]]]*b[[op[2L]]],tolerance=2e-13)
            expect_equal(unname(z$p.kw[,,1L]),a$normal*b[[op[2L]]],tolerance=2e-13)
            expect_equal(unname(z$p.kw[,,2L]),a[[op[1L]]]*b$normal,tolerance=2e-13)
          }
        }
      }
    })
    test_that("positive tied radii and extended NN retain dense/tree agreement", {
      train <- data.frame(x=seq_len(40L)/4)
      at <- data.frame(x=c(.25,2.25,4.5,8))
      for(type in c("generalized_nn","adaptive_nn")) for(k in c(6L,80L)) {
        values <- lapply(list(FALSE,TRUE,"auto"),function(tree) {
          options(np.tree=tree)
          npksum(txdat=train,exdat=at,bws=k,bwtype=type,ckertype="uniform",
                 return.kernel.weights=TRUE)$kw
        })
        expect_identical(values[[1L]],values[[2L]])
        expect_identical(values[[1L]],values[[3L]])
      }
    })
  }
  if(exists(".npRmpi_with_local_regression",mode="function"))
    .npRmpi_with_local_regression(run()) else run()
})
