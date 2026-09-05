test_that("least-squares quantile generics preserve formula subset evaluation", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  dat <- data.frame(x = seq(.05, .95, length.out = 24),
                    y = sin(seq_len(24) / 3) + seq_len(24) / 24)
  dat$keep <- dat$x > .25
  fo <- as.formula("y ~ x", env = baseenv())
  for (family in c("nplsqregbw", "nplsqreg")) {
    fun <- get(family, mode = "function")
    direct <- get(paste0(family, ".formula"), environment(fun), inherits = TRUE)
    for (tau in list(.5, c(.25, .75))) {
      selected <- dat[dat$keep, , drop = FALSE]
      controls <- list(scale = rep(1, nrow(selected)), tau = tau,
        bandwidth.compute = FALSE, nomad = FALSE, regtype = "ll",
        nmulti = 1L, itmax = 5L, random.seed = 17L)
      reference <- do.call(fun, c(if (family == "nplsqregbw")
        list(xdat = selected["x"], ydat = selected$y) else
        list(txdat = selected["x"], tydat = selected$y), controls))
      payload <- function(z) {
        b <- if (family == "nplsqregbw") z else z$bws
        list(xdat = b$xdat, ydat = b$ydat, bandwidth = b$reg.bws$bw,
             fitted = if (family == "nplsqreg") fitted(z) else NULL)
      }
      keep <- !dat$keep
      for (expr in list(quote(x > .25), quote(keep))) {
        actual <- eval(as.call(c(list(quote(fun)),
          list(bws = quote(fo), data = quote(dat), subset = expr), controls)))
        expect_identical(payload(actual), payload(reference))
        direct.result <- eval(as.call(c(list(quote(direct)),
          list(bws = quote(fo), data = quote(dat), subset = expr), controls)))
        expect_identical(payload(actual), payload(direct.result))
      }
      data.calls <- 0L
      counter <- new.env(parent = emptyenv())
      counter$calls <- 0L
      counter.fo <- as.formula("y ~ x", env = list2env(list(counter = counter), parent = baseenv()))
      actual <- eval(as.call(c(list(quote(fun)), list(bws = quote(counter.fo),
        data = quote({data.calls <- data.calls + 1L; dat}),
        subset = quote({counter$calls <- counter$calls + 1L; x > .25})), controls)))
      expect_identical(data.calls, 1L)
      expect_identical(counter$calls, 1L)
      expect_identical(payload(actual), payload(reference))
      partial <- eval(as.call(c(list(quote(fun)),
        list(bws = quote(fo), data = quote(dat), sub = quote(x > .25)), controls)))
      expect_identical(payload(partial), payload(reference))
    }
  }
})

test_that("least-squares formula subsets retain omission and environment policy", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  dat <- data.frame(x = seq(.05, .95, length.out = 24), y = sin(seq_len(24)))
  fo <- as.formula("y ~ x", env = baseenv())
  selections <- list(NULL, seq_len(24) > 3, 4:24, c(NA_integer_, 4:24))
  for (family in c("nplsqregbw", "nplsqreg")) {
    fun <- get(family, mode = "function")
    direct <- get(paste0(family, ".formula"), environment(fun), inherits = TRUE)
    for (selection in selections) {
      selected <- do.call(model.frame, list(formula = fo, data = dat,
        subset = selection, na.action = na.omit))
      controls <- list(scale = rep(1, nrow(selected)), bandwidth.compute = FALSE,
        nomad = FALSE, regtype = "ll", nmulti = 1L, itmax = 5L, random.seed = 17L)
      args <- c(list(bws = fo, data = dat, subset = selection, na.action = na.omit), controls)
      actual <- do.call(fun, args)
      reference <- do.call(direct, args)
      if (family == "nplsqregbw") {
        expect_identical(actual$xdat, reference$xdat)
        expect_identical(actual$ydat, reference$ydat)
        expect_identical(actual$reg.bws$bw, reference$reg.bws$bw)
      } else {
        expect_identical(fitted(actual), fitted(reference))
        expect_identical(predict(actual, newdata = data.frame(x = c(.2, .5))),
                         predict(reference, newdata = data.frame(x = c(.2, .5))))
      }
    }
    local.keep <- rep(TRUE, nrow(dat))
    expect_error(fun(fo, data = dat, subset = local.keep,
      scale = rep(1, nrow(dat)), bandwidth.compute = FALSE), "local.keep.*not found")
    expect_error(fun(fo, data = dat, subset = stop("subset sentinel"),
      scale = rep(1, nrow(dat)), bandwidth.compute = FALSE), "subset sentinel")
  }
})

test_that("least-squares dispatch retains non-subset dots and promise caching", {
  calls <- 0L
  forward <- function(...) {
    args <- .nplsqreg_formula_dispatch_args(...)
    all <- list(...)
    list(args = args, all = all)
  }
  z <- forward(txdat = {calls <- calls + 1L; NULL}, subset = TRUE,
               unnamed = 3, FALSE)
  expect_identical(calls, 1L)
  expect_identical(z$args, z$all[-2L])
  expect_identical(.nplsqreg_formula_dispatch_args(data = 1,
    subset = stop("must remain lazy")), list(data = 1))
  expect_identical(.nplsqreg_formula_dispatch_args(1, sub = stop("lazy")), setNames(list(1), ""))
})
