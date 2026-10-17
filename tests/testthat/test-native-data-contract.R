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
