# Heavy MPI delete-one oracle; each explicit selector runs in its own bounded pool.
npRmpi_test_gnn_deleted_query_count <- function(trees, kernels) {
  old <- options(np.messages=FALSE, np.extendednn=TRUE)
  on.exit(options(old), add=TRUE)
  x <- data.frame(x=c(.03,.09,.18,.29,.43,.57,.71,.84,.97))
  y <- .15+.55*x$x+.1*sin(14*x$x)
  yd <- data.frame(y=y); n <- nrow(x)
  for (tree in trees) for (kernel in kernels) {
    options(np.tree=tree)
    bounds <- if (kernel=="beta") list(ckerbound="fixed", ckerlb=0, ckerub=1) else list()
    cbounds <- if (kernel=="beta") list(cxkerbound="fixed", cxkerlb=0, cxkerub=1,
      cykerbound="fixed", cykerlb=0, cykerub=1) else list()
    for (k in c(n-3L, n-2L, n-1L, n+2L)) {
      info <- paste(tree, kernel, k)
      ub <- do.call(npudensbw, c(list(dat=x, bws=k, bwtype="generalized_nn",
        bwmethod="cv.ml", ckertype=kernel, bandwidth.compute=FALSE), bounds))
      ud <- vapply(seq_len(n), function(i) fitted(npudens(bws=ub,
        tdat=x[-i,,drop=FALSE], edat=x[i,,drop=FALSE])), 0.)
      expect_equal(as.numeric(npudensbw.bandwidth(dat=x, bws=ub, eval.only=TRUE,
        nmulti=1L, invalid.penalty="dbmax")$fval), sum(log(ud)),
        tolerance=2e-12, info=info)
      for (rt in c("lc", "ll", "lp")) {
        rb <- do.call(npregbw, c(list(xdat=x, ydat=y, bws=k,
          bwtype="generalized_nn", regtype=rt, ckertype=kernel,
          bandwidth.compute=FALSE), if(rt=="lp") list(degree=2L) else list(), bounds))
        pred <- vapply(seq_len(n), function(i) fitted(npreg(bws=rb,
          txdat=x[-i,,drop=FALSE], tydat=y[-i], exdat=x[i,,drop=FALSE])), 0.)
        expect_equal(as.numeric(.npregbw_eval_only(x,y,rb,
          invalid.penalty="dbmax")$objective), mean((y-pred)^2),
          tolerance=2e-12, info=paste(info,rt))
        # Generic beta/extended MPI hats require the separate H adapter.
        # Keep their literal oracle in the campaign dependency probe.
        if (kernel != "beta" && k <= n-1L)
          expect_equal(as.numeric(npreghat(rb, txdat=x, y=y, output="apply",
          leave.one.out=TRUE)), pred, tolerance=2e-12, info=paste(info,rt))
      }
      cb <- do.call(npcdensbw, c(list(xdat=x, ydat=yd, bws=c(k,k),
        bwtype="generalized_nn", bwmethod="cv.ml", regtype="lc",
        cxkertype=kernel, cykertype=kernel, bandwidth.compute=FALSE), cbounds))
      cd <- vapply(seq_len(n), function(i) fitted(npcdens(bws=cb,
        txdat=x[-i,,drop=FALSE], tydat=yd[-i,,drop=FALSE],
        exdat=x[i,,drop=FALSE], eydat=yd[i,,drop=FALSE])), 0.)
      expect_equal(as.numeric(.npcdensbw_eval_only(x,yd,cb,
        invalid.penalty="dbmax")$objective), sum(log(cd)),
        tolerance=2e-12, info=info)
      if (k==n-1L) {
        options(np.extendednn=FALSE)
        expect_identical(as.numeric(npudensbw.bandwidth(dat=x, bws=ub,
          eval.only=TRUE, nmulti=1L, invalid.penalty="dbmax")$fval),
          -.Machine$double.xmax)
        expect_identical(as.numeric(.npregbw_eval_only(x,y,rb,
          invalid.penalty="dbmax")$objective), .Machine$double.xmax)
        expect_identical(as.numeric(.npcdensbw_eval_only(x,yd,cb,
          invalid.penalty="dbmax")$objective), -.Machine$double.xmax)
        if (kernel != "beta")
          expect_error(npreghat(rb,txdat=x,y=y,output="apply",leave.one.out=TRUE))
        options(np.extendednn=TRUE)
      }
    }
  }
}
