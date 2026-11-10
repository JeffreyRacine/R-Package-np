test_that("native frame argument validation inspects names without forcing promises", {
  reject <- getFromNamespace(".np_reject_native_frame_args", 'npRmpi')
  for (arg in c("subset", "na.action")) {
    expressions <- as.list(substitute(list(stop("must not evaluate"))))[-1L]
    names(expressions) <- arg
    expect_error(reject(expressions), "supported only with the formula interface", fixed=TRUE)
  }
  expect_null(reject(alist(ckertype="gaussian", bandwidth.compute=FALSE)))
})

test_that("every native constructor leaf rejects ignored frame arguments", {
  leafs <- c("npregbw.NULL", "npregbw.default", "npregbw.rbandwidth",
             "npudensbw.NULL", "npudensbw.default", "npudensbw.bandwidth",
             "npudistbw.NULL", "npudistbw.default", "npudistbw.dbandwidth",
             "npindexbw.NULL", "npindexbw.default", "npindexbw.sibandwidth",
             "npplregbw.NULL", "npplregbw.default", "npplregbw.plbandwidth",
             "npscoefbw.NULL", "npscoefbw.default", "npscoefbw.scbandwidth",
             "npksum.integer", "npksum.numeric", "npksum.default")
  for (leaf in leafs) for (arg in c("subset", "na.action")) {
    fn <- getFromNamespace(leaf, 'npRmpi')
    supplied <- alist(xdat=1:12, dat=1:12, txdat=1:12, ydat=seq_len(12),
                     tydat=seq_len(12), zdat=1:12, bws=.3)
    supplied <- supplied[names(supplied) %in% names(formals(fn))]
    supplied[[arg]] <- quote(stop("must not evaluate frame argument"))
    call <- as.call(c(list(fn), supplied))
    expect_error(eval(call), "supported only with the formula interface", fixed=TRUE, info=paste(leaf,arg))
    supplied[[arg]] <- NULL
    # Explicit NULL is still an unsupported native argument.
    supplied[arg] <- list(NULL)
    expect_error(eval(as.call(c(list(fn),supplied))), "supported only with the formula interface", fixed=TRUE, info=paste(leaf,arg,"NULL"))
  }
})
