test_that("recoding distinguishes NA category labels from missing codes", {
  lev <- c("a", "b", NA_character_)
  for (ordered in c(FALSE, TRUE)) {
    f <- factor(c("a", NA, "b", "a"), levels = lev, exclude = NULL,
                ordered = ordered)
    is.na(f)[4] <- TRUE
    for (out in list(.np_factor_with_levels(f, lev), cast(f, f),
                     adjustLevels(data.frame(f = f), untangle(data.frame(f = f)))$f)) {
      expect_identical(as.integer(out), c(1L, 3L, 2L, NA_integer_))
      expect_identical(levels(out), lev)
      expect_identical(is.ordered(out), ordered)
    }
    plain <- factor(c("a","b","a"), ordered = ordered)
    expect_identical(.np_factor_with_levels(plain, levels(plain)),
                     factor(plain, levels = levels(plain)))
    grid <- .np_plot_conmode_grid_values(f, 3, c(0,1))
    expect_identical(as.integer(grid), 1:3)
    expect_identical(as.integer(.np_plot_conmode_cast_like(f, f)),
                     as.integer(f))
    endpoints <- npCategoricalFirstDifferenceFrames(data.frame(f = f), 1, "test")
    expect_identical(is.na(endpoints$lower$f), is.na(f))
    expect_identical(is.na(endpoints$upper$f), is.na(f))
    expect_false(is.na(endpoints$upper$f[2]))
    selected <- .npConmodeSelectedFactor(c(1L,3L,0L), lev, ordered)
    expect_identical(as.integer(selected), c(1L,3L,NA_integer_))
    padded <- .npConmodeNapredictRows(2L, selected)
    expect_identical(as.integer(padded), c(1L,NA_integer_,3L,NA_integer_))
  }
})

test_that("NA-labelled categories satisfy direct weights and relabelling oracles", {
  old <- options(np.messages=FALSE,np.largeh=FALSE,np.largelambda=FALSE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(x=factor(rep(c("a","b",NA),20),exclude=NULL))
  y <- sin(seq_len(nrow(x)))
  lambda <- .1
  code <- as.integer(x$x)
  K <- outer(code,code,function(a,b) ifelse(a==b,1-lambda,lambda/2))
  expect_equal(as.double(fitted(npudens(tdat=x,bws=lambda))),rowMeans(K),
               tolerance=1e-13)
  expect_equal(as.double(fitted(npreg(txdat=x,tydat=y,bws=lambda))),
               as.vector(K%*%y/rowSums(K)),tolerance=1e-13)
  for (ordered in c(FALSE,TRUE)) {
    xx <- x
    if (ordered) class(xx$x) <- c("ordered","factor")
    renamed <- xx
    levels(renamed$x) <- c("a","b","missing-label")
    for (family in c("npreg","npudens","npcdens","npcdist")) {
      args <- if (family=="npudens") list(tdat=xx,bws=lambda,se=TRUE) else
        list(txdat=xx,tydat=y,bws=if(family=="npreg") lambda else c(.4,lambda),
             se=TRUE,gradients=TRUE)
      actual <- do.call(get(family),args)
      args[[if(family=="npudens") "tdat" else "txdat"]] <- renamed
      oracle <- do.call(get(family),args)
      expect_equal(fitted(actual),fitted(oracle),tolerance=1e-12)
      expect_equal(se(actual),se(oracle),tolerance=1e-12)
      if(family!="npudens")
        expect_equal(gradients(actual),gradients(oracle),tolerance=1e-12)
    }
    support <- .np_entropy_factor_support(xx$x)
    expect_identical(as.integer(support$x),as.integer(xx$x))
    expect_equal(sort(as.integer(support$evaluation)),1:3)
  }
})
