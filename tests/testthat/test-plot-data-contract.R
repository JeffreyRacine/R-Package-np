test_that("plot routes reject data they do not consume", {
  msg <- "This plot method does not consume `data=` on the selected route. Refit the model with the intended data and plot that model, or supply the supported explicit plotting inputs."
  engines <- c("bandwidth", "dbandwidth", "rbandwidth", "conbandwidth",
               "condbandwidth", "scbandwidth", "sibandwidth")
  for (kind in engines) {
    owner <- get(paste0(".np_plot_", kind, "_engine"))
    for (call in list(quote(owner(data = NULL)),
                      quote(owner(data = stop("data was forced"))),
                      quote(owner(data = )))) {
      expect_identical(tryCatch(eval(call), error = conditionMessage), msg)
    }
    for (expr in list(quote(list(data = NULL)),
                      quote(list(data = stop("data was forced"))),
                      quote(list(data = )))) {
      expect_identical(tryCatch(.np_plot_validate_public_dots(
        expr[-1L], method = owner), error = conditionMessage), msg)
    }
  }
  for (call in list(quote(.np_plot_plbandwidth_engine(bws = list(), data = NULL)),
                    quote(.np_plot_plbandwidth_engine(bws = list(), data = stop("data was forced"))),
                    quote(.np_plot_plbandwidth_engine(bws = list(), data = )))) {
    expect_identical(tryCatch(eval(call), error = conditionMessage), msg)
  }
  for (call in list(quote(.np_plot_plbandwidth_engine(bws = list(formula = y ~ x | z), xdat = 1, ydat = 1, zdat = 1, data = NULL)),
                    quote(.np_plot_plbandwidth_engine(bws = list(formula = y ~ x | z), xdat = 1, ydat = 1, zdat = 1, data = stop("data was forced"))),
                    quote(.np_plot_plbandwidth_engine(bws = list(formula = y ~ x | z), xdat = 1, ydat = 1, zdat = 1, data = )))) {
    expect_identical(tryCatch(eval(call), error = conditionMessage), msg)
  }
})

test_that("explicit data-formal suppliedness respects R matching", {
  f <- function(x, data = NULL, ...) .np_data_formal_supplied(sys.call(), sys.function())
  expect_false(f())
  expect_false(f(newdata = stop("not evaluated")))
  expect_true(f(data = NULL))
  expect_true(f(data = ))
  expect_true(f(dat = ))
  expect_true(f(1, ))
  expect_true(f(data = stop("not evaluated")))
  forward <- function(...) f(...)
  nested <- function(...) forward(...)
  expect_false(nested(1))
  expect_false(nested(1, newdata = stop("not evaluated")))
  expect_true(nested(1, data = NULL))
  expect_true(nested(1, data = ))
  expect_true(nested(1, dat = ))
  expect_true(nested(1, ))
  expect_true(nested(1, data = stop("not evaluated")))
  expect_true(do.call(nested, list(1, data = NULL)))
})
