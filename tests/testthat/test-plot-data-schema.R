test_that("plot data has explicit coordinates and absolute endpoints", {
  grid <- data.frame(lower = ordered(c("low", "high"), levels = c("low", "high")))
  values <- npRmpi:::.np_plot_data_values(c(1, 2), "mean", c(1.25, 2.5), c(.1, .2))
  # The corrected estimate need not be the interval midpoint.
  intervals <- npRmpi:::.np_plot_data_intervals(
    list(simultaneous = cbind(c(.8, 1.5), c(1.4, 2.8))), "mean", .05, "f1")
  panel <- npRmpi:::.np_plot_data_panel(grid, "X", values, intervals)
  targets <- npRmpi:::.np_plot_data_target("mean", "mean", inference_id = "i1")
  families <- list(f1 = list(event = "all target rows covered", scope = "finite grid",
    members = data.frame(panel_id = "p1", row_id = 1:2, target_id = "mean")))
  out <- npRmpi:::.np_plot_data("npreg", list(), list(p1 = panel), targets,
    list(i1 = list(errors = "bootstrap")), families)
  expect_identical(class(out), "np_plot_data")
  expect_identical(out$panels$p1$grid$eval_1, grid$lower)
  expect_identical(out$panels$p1$variables$name, "lower")
  # Resolve directly here until public registration is qualified.
  flat <- npRmpi:::as.data.frame.np_plot_data(out)
  expect_identical(flat$lower, c(.8, 1.5))
  expect_identical(flat$upper, c(1.4, 2.8))
  expect_identical(flat$bias_corrected, c(1.25, 2.5))
  expect_false(any(c("merr", "inference.bands", "bws") %in% names(out)))
  expect_error(npRmpi:::as.data.frame.np_plot_data(out, panel = 0), "invalid panel")
  bad <- out; bad$panels$p1$intervals$lower[1] <- 9
  expect_error(npRmpi:::.np_plot_data_validate(bad), "reversed interval")
  bad <- out; bad$panels$p1$intervals$row_id[2] <- 1L
  expect_error(npRmpi:::.np_plot_data_validate(bad), "duplicate keys")
  bad <- out; bad$panels$p1$intervals$target_id[1] <- "wrong"
  expect_error(npRmpi:::.np_plot_data_validate(bad), "invalid row/target")
})

test_that("singleton axes and unavailable rows remain explicit", {
  values <- npRmpi:::.np_plot_data_values(c(NA_real_), "t")
  intervals <- npRmpi:::.np_plot_data_empty_intervals()
  panel <- npRmpi:::.np_plot_data_panel(data.frame(x = 3), "X", values, intervals)
  out <- npRmpi:::.np_plot_data("npudens", list(), list(p1 = panel),
    npRmpi:::.np_plot_data_target("t", "density"))
  flat <- npRmpi:::as.data.frame.np_plot_data(out)
  expect_equal(nrow(flat), 1L)
  expect_identical(flat$row_id, 1L)
  expect_identical(flat$status, "undefined")
  expect_false("lower" %in% names(flat))
  expect_error(npRmpi:::as.data.frame.np_plot_data(out, band = "pointwise"), "no intervals")
  out$panels$p2 <- panel
  expect_error(npRmpi:::as.data.frame.np_plot_data(out), "select panel explicitly")
})

test_that("endpoint orientation is explicit rather than guessed from signs", {
  reference <- c(2, -1)
  offsets <- cbind(c(.1, -.3), c(.4, .2))
  widths <- cbind(c(-.1, .3), c(.4, .2))
  endpoints <- cbind(c(2.1, -1.3), c(2.4, -.8))
  expect_equal(npRmpi:::.np_plot_data_endpoints(offsets, reference, "offsets"), endpoints)
  expect_equal(npRmpi:::.np_plot_data_endpoints(widths, reference, "widths"), endpoints)
  expect_identical(npRmpi:::.np_plot_data_endpoints(endpoints, reference, "endpoints"), endpoints)
  expect_false(isTRUE(all.equal(
    npRmpi:::.np_plot_data_endpoints(offsets, reference, "widths"), endpoints)))
})

# These pure tests cover schema and coordinate identity only; numerical method
# coverage remains in the existing family tests rather than a new heavy lane.
test_that("plot metadata preserves contrast direction and coordinate roles", {
  ord <- data.frame(x = ordered(c("low", "middle", "high"),
                                levels = c("low", "middle", "high")))
  z <- npRmpi:::.np_plot_data_contrast(ord, 1L)
  expect_identical(z$lower_level, c("low", "low", "middle"))
  expect_identical(z$upper_level, c("middle", "middle", "high"))
  expect_identical(z$operation, "upper minus lower")
  ord$x <- factor(ord$x, levels = c("middle", "low", "high"), ordered = FALSE)
  z <- npRmpi:::.np_plot_data_contrast(ord, 1L)
  expect_identical(z$lower_level, rep("middle", 3L))
  expect_identical(z$upper_level, c("low", "middle", "high"))
  # Repaired display names must not confuse X and Y or two X coordinates.
  observed <- data.frame(a = 1:3, b = 4:6, c = 7:9)
  names(observed) <- c("x.one", "x.one.1", "x.one.2")
  grid <- observed; names(grid) <- c("x one", "x.one", "x one")
  grid[[1L]] <- 0L
  raw <- npRmpi:::.np_plot_data_publish_display(list(), observed, c("X", "X", "Y"),
                                            response = 11:13, overlay = TRUE, rug = TRUE)
  z <- npRmpi:::.np_plot_data_display(raw, list(grid = grid, roles = c("X", "X", "Y")), list())
  expect_identical(z$overlay$grid$eval_2, 4:6)
  expect_identical(z$overlay$grid$eval_3, 7:9)
  expect_identical(z$rug, z$overlay$grid)
})

test_that("plot schema rejects field and inference ownership corruption", {
  panel <- npRmpi:::.np_plot_data_panel(data.frame(x = 1), "X",
    npRmpi:::.np_plot_data_values(2, "t"), npRmpi:::.np_plot_data_empty_intervals())
  target <- npRmpi:::.np_plot_data_target("t", "mean", inference_id = "i")
  out <- npRmpi:::.np_plot_data("npreg", list(), list(p = panel), target, list(i = list()))
  bad <- out; bad$targets$inference_id <- "missing"
  expect_error(npRmpi:::.np_plot_data_validate(bad), "undefined inference record")
  bad <- out; bad$panels$p$values$estimate <- "2"
  expect_error(npRmpi:::.np_plot_data_validate(bad), "invalid fields")
  bad <- out; bad$obsolete <- list()
  expect_error(npRmpi:::.np_plot_data_validate(bad), "unexpected top-level")
})


test_that("coefficient inference admission agrees for fitted and bandwidth objects", {
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves())
  on.exit(close_mpi_slaves(force=TRUE), add=TRUE)
  old <- options(np.messages=FALSE)
  on.exit(options(old), add=TRUE)
  # Keep fitted call environments separate from the changing test frame, so
  # full serialization checks measure model mutation, not test bookkeeping.
  fixture <- evalq({
    set.seed(94114)
    x <- data.frame(x=runif(24L)); z <- data.frame(z=runif(24L))
    y <- 1.2*x$x + cos(2*pi*z$z) + rnorm(24L, sd=.05)
    sb <- npscoefbw(xdat=x, zdat=z, ydat=y, bws=.4,
                    bandwidth.compute=FALSE, regtype="ll")
    sf <- npscoef(bws=sb, txdat=x, tzdat=z, tydat=y, betas=TRUE)
    pb <- npplregbw(xdat=x, zdat=z, ydat=y, bws=matrix(.4,2L,1L),
                    bandwidth.compute=FALSE)
    pf <- npplreg(bws=pb)
    list(x=x,z=z,y=y,cases=list(list(bw=sb,fit=sf),list(bw=pb,fit=pf)))
  }, new.env(parent=globalenv()))
  x <- fixture$x; z <- fixture$z; y <- fixture$y
  for (case in fixture$cases) {
    saved <- serialize(case$fit,NULL); rng <- .Random.seed
    for (object in list(case$bw,case$fit)) {
      common <- list(x=object,xdat=x,zdat=z,ydat=y,coef=TRUE,
                     output="data",perspective=FALSE,neval=5L)
      p <- do.call(plot,c(common,list(errors="none")))
      expect_s3_class(p,"np_plot_data")
      expect_true(all(p$targets$quantity=="coefficient"))
      expect_error(do.call(plot,c(common,list(errors="asymptotic",band="simultaneous"))),
                   "coefficient", fixed=TRUE)
      expect_error(do.call(plot,c(common,list(errors="bootstrap",B=12L))),
                   "coefficient", fixed=TRUE)
      if(inherits(case$bw,"scbandwidth"))
        expect_error(do.call(plot,c(common,list(errors="bootstrap"))),
                     "smooth coefficient displays", fixed=TRUE)
    }
    expect_identical(saved,serialize(case$fit,NULL))
    expect_identical(rng,.Random.seed)
  }
})
