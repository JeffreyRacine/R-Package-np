# The policy and expected diagnostic are specified here independently of the
# implementation helper. These checks do not enter numerical/MPI work.
test_that("selected native owners reject data without forcing its value", {
  owners <- c(
    "npcdensbw.conbandwidth",
    "npcdensbw.NULL",
    "npcdensbw.default",
    "npcdens.call",
    "npcdens.conbandwidth",
    "npcdistbw.condbandwidth",
    "npcdistbw.NULL",
    "npcdistbw.default",
    "npcdist.call",
    "npcdist.condbandwidth",
    "npconmode.call",
    "npconmode.conbandwidth",
    "npudensbw.NULL",
    "npudensbw.bandwidth",
    "npudensbw.default",
    "npudens.call",
    "npudens.bandwidth",
    "npudistbw.NULL",
    "npudistbw.dbandwidth",
    "npudistbw.default",
    "npudist.call",
    "npudist.dbandwidth",
    "npksum.integer",
    "npksum.numeric",
    "npksum.default",
    "nplsqregbw.lsqregressionbandwidth",
    "nplsqregbw.NULL",
    "nplsqregbw.default",
    "nplsqreg.lsqregressionbandwidth",
    "nplsqreg.NULL",
    "nplsqreg.default",
    "npplregbw.NULL",
    "npplregbw.plbandwidth",
    "npplregbw.default",
    "npplreg.call",
    "npplreg.plbandwidth",
    "npqreg.call",
    "npqreg.condbandwidth",
    "npreghat.call",
    "npreghat.npregression",
    "npreghat.rbandwidth",
    "npregbw.NULL",
    "npregbw.rbandwidth",
    "npregbw.default",
    "npreg.call",
    "npreg.rbandwidth",
    "npsigtest.call",
    "npsigtest.rbandwidth",
    "npindexbw.NULL",
    "npindexbw.default",
    "npindexbw.sibandwidth",
    "npindex.call",
    "npindex.sibandwidth",
    "npscoefbw.NULL",
    "npscoefbw.scbandwidth",
    "npscoefbw.default",
    "npscoef.call",
    "npscoef.scbandwidth",
    "npregiv.default",
    "npregivderiv.default")
  diagnostic <- "`data=` is supported only with the formula interface. For native calls, supply vectors, matrices or data frames directly, or use `with(data, ...)`."
  for (owner in owners) {
    fun <- get(owner, envir=environment(npreg), inherits=FALSE)
    for (expr in list(quote(fun(data=NULL)),
                      quote(fun(data=stop("data promise was forced"))),
                      quote(fun(data=)))) {
      actual <- tryCatch(eval(expr), error=conditionMessage)
      expect_identical(actual, diagnostic, info=owner)
    }
  }
})

test_that("IV native bare-name misuse has an actionable data diagnostic", {
  for (fun in list(npregiv,npregivderiv)) {
    expect_error(fun(y=.np_absent_iv_y,z=.np_absent_iv_z,w=.np_absent_iv_w,
                     data=stop("data promise was forced")),
                 "`data=` is supported only with the formula interface",fixed=TRUE)
  }
})

test_that("CMS native data suppliedness includes NULL and partial formal names", {
  for(fun in list(npcmstest,npqcmstest)) {
    expect_error(fun(xdat=1:8,ydat=1:8,data=NULL),
                 "`data=` is supported only with the formula interface",fixed=TRUE)
    expect_error(fun(xdat=1:8,ydat=1:8,dat=stop("data promise was forced")),
                 "`data=` is supported only with the formula interface",fixed=TRUE)
  }
})

# This closes forwarded formal matching, distinct from the dots-name guards.
# These public calls must stop before model validation, fitting or MPI dispatch.
test_that("CMS data rejection survives wrapper forwarding", {
  forward <- function(fun, ...) fun(...)
  nested <- function(fun, ...) forward(fun, ...)
  diagnostic <- "`data=` is supported only with the formula interface"
  for (fun in list(npcmstest, npqcmstest)) {
    expect_error(nested(fun, xdat = 1:8, ydat = 1:8, data = NULL),
                 diagnostic, fixed = TRUE)
    expect_error(nested(fun, xdat = 1:8, ydat = 1:8, data = ),
                 diagnostic, fixed = TRUE)
    expect_error(nested(fun, xdat = 1:8, ydat = 1:8,
                       dat = stop("data promise was forced")),
                 diagnostic, fixed = TRUE)
  }
})

# Public first-force coverage complements the selected-owner tests above: those
# owners are never reached if an unbound observation name fails during dispatch.
test_that("public bare native names receive data-interface guidance", {
  roles <- list(
    npregbw = c("xdat", "ydat"), npudensbw = "dat", npudistbw = "dat",
    npcdensbw = c("xdat", "ydat"), npcdistbw = c("xdat", "ydat"),
    npindexbw = c("xdat", "ydat"), npscoefbw = c("xdat", "ydat", "zdat"),
    npplregbw = c("xdat", "ydat", "zdat"), nplsqregbw = c("xdat", "ydat"),
    npreg = c("tydat", "txdat"), npudens = "tdat", npudist = "tdat",
    npcdens = c("txdat", "tydat"), npcdist = c("txdat", "tydat"),
    npindex = c("txdat", "tydat"), npscoef = c("txdat", "tydat", "tzdat"),
    npplreg = c("txdat", "tydat", "tzdat"), nplsqreg = c("txdat", "tydat"),
    npqreg = c("txdat", "tydat"), npconmode = c("txdat", "tydat"),
    npksum = "txdat", npreghat = c("txdat", "y"),
    npsigtest = c("xdat", "ydat"), npcmstest = c("xdat", "ydat"),
    npqcmstest = c("xdat", "ydat"), npregiv = c("y", "z", "w"),
    npregivderiv = c("y", "z", "w"))
  for (family in names(roles)) {
    for (reversed in c(FALSE, TRUE)) {
      names <- if (reversed) rev(roles[[family]]) else roles[[family]]
      inputs <- setNames(lapply(names, function(x) as.name(paste0(".np_unbound_", x))), names)
      inputs["data"] <- list(quote(NULL))
      fun <- get(family, envir = environment(npreg), inherits = FALSE)
      expr <- as.call(c(list(fun), inputs))
      expect_error(eval(expr), "`data=` is supported only with the formula interface",
                   fixed = TRUE, info = paste(family, reversed))
    }
  }
})

test_that("early data diagnostics preserve factories and active bindings", {
  factory <- function() stop("original user factory error", call. = FALSE)
  expect_error(npreg(txdat = factory(), tydat = 1:4, data = NULL),
               "^original user factory error$")
  scope <- new.env(parent = environment())
  makeActiveBinding("active", factory, scope)
  expect_error(eval(quote(npreg(txdat = active, tydat = 1:4, data = NULL)), scope),
               "^original user factory error$")
  # A formula itself must still resolve through ordinary R evaluation.
  expect_error(npreg(formula = .np_unbound_formula, data = NULL),
               "^object '.np_unbound_formula' not found$")
})

# R retains a forced promise even if a previous argument removes its source name.
test_that("early diagnostic preserves an already cached observation promise", {
  scope <- new.env(parent = environment())
  scope$x <- 1:4
  materialize <- function(...) {
    .np_formula_dispatch_args(NULL, substitute(list(...))[-1L], environment(),
      data.roles = "txdat", call = sys.call(), caller = parent.frame())
  }
  scope$materialize <- materialize
  scope$cache <- function() {
    for (i in seq_len(sys.nframe())) if (identical(sys.function(i), materialize)) {
      eval(quote(..2), envir = sys.frame(i))
      break
    }
    rm("x", envir = scope)
    NULL
  }
  result <- eval(quote(materialize(data = cache(), txdat = x)), scope)
  expect_identical(result, list(data = NULL, txdat = 1:4))
  expect_false(exists("x", envir = scope, inherits = FALSE))
})
