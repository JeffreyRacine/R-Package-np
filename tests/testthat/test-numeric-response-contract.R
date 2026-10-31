test_that("numeric response validation is non-coercing and shape-specific", {
  good <- list(1:4, c(1, NA_real_, NaN, Inf), setNames(1:4, letters[1:4]),
               numeric(), ts(1:4))
  for (y in good) {
    before <- y
    expect_null(.np_require_numeric_response(y, "probe", "y"))
    expect_identical(y, before)
  }
  bad <- list(factor(1:4), ordered(c(10, 30, 50)), c("1", "2"), TRUE,
              1+1i, as.raw(1:4), list(1, 2), NULL, matrix(1:4, 4),
              array(1:4, c(2, 2, 1)), data.frame(y=1:4), Sys.Date())
  for (y in bad)
    expect_error(.np_require_numeric_response(y, "probe", "y"),
                 "must be a numeric \\(integer or double\\) response vector")
})

test_that("semiparametric public outcomes are numeric and predictors keep their roles", {
  withr::local_options(np.messages = FALSE)

  n <- 24L
  d <- data.frame(x = seq(-1, 1, length.out=n), w = sin(seq_len(n)),
                  z = cos(seq_len(n)/3), y = sin(seq_len(n)/4))
  f <- factor(rep(c("50", "10", "30"), 8), levels=c("50", "10", "30", "70"))
  o <- ordered(rep(c("30", "10", "50"), 8), levels=c("10", "30", "50", "70"))
  constructors <- list(
    plreg = function(y) npplregbw(xdat=d[c("x","w")], ydat=y, zdat=d["z"],
                      bws=matrix(.6, 3, 1), bandwidth.compute=FALSE),
    scoef = function(y) npscoefbw(xdat=d[c("x","w")], ydat=y, zdat=d["z"],
                      bws=.6, bandwidth.compute=FALSE),
    index = function(y) npindexbw(xdat=d[c("x","w")], ydat=y,
                      bws=c(1,.2,.6), bandwidth.compute=FALSE))
  fitters <- list(plreg=npplreg, scoef=npscoef, index=npindex)
  for (family in names(constructors)) {
    b <- constructors[[family]](d$y)
    for (bad in list(f,o)) {
      expect_error(constructors[[family]](bad), "must be a numeric")
      args <- list(bws=b, txdat=d[c("x","w")], tydat=bad, se=FALSE)
      if (family != "index") args$tzdat <- d["z"]
      expect_error(do.call(fitters[[family]],args), "must be a numeric")
      args$tydat <- d$y; args$eydat <- bad
      expect_error(do.call(fitters[[family]],args), "must be a numeric")
    }
    fit <- fitters[[family]](b, se=FALSE)
    expect_length(fitted(fit), n)
  }
  # These are predictor responses owned by npplreg, never primary outcomes.
  for (cat in list(f,o)) {
    b <- npplregbw(xdat=data.frame(cat=cat), ydat=d$y, zdat=d["z"],
                   bws=matrix(.6,2,1), bandwidth.compute=FALSE)
    expect_identical(b$bw[[2L]][[".np.native.training"]]$ydat, cat)
    expect_length(fitted(npplreg(b)), n)
    expect_length(fitted(npreg(unserialize(serialize(b$bw[[2L]],NULL)))), n)
  }
})

test_that("optional Ichimura scores do not turn unused outcomes into errors", {
  y <- factor(c("1", "2"))
  expect_null(.np_index_diagnostics_response(y, "ichimura"))
  expect_error(.np_index_diagnostics_response(y, "ichimura", required=TRUE),
               "must be a numeric")
  expect_error(.np_index_diagnostics_response(NULL, "ichimura", required=TRUE),
               "must be a numeric")
  expect_identical(.np_index_diagnostics_response(1:2, "ichimura"), 1:2)
  expect_error(.np_index_diagnostics_response(y, "kleinspady", required=TRUE),
               "Klein and Spady")
})
