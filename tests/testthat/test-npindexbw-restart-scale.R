test_that("automatic optim restarts use the drawn beta and preserve the U draw", {
  controls <- .npindexbw_h_start_controls(0.7, 0.9, 1.2, 0)
  for (p in c(1L, 2L, 3L, 5L)) {
    i <- seq_len(48L)
    x <- vapply(seq_len(p), function(j) sin(i * sqrt(j + 1)) * j, numeric(48L))
    beta <- if (p > 1L) seq_len(p - 1L) / 3 else numeric()
    first <- as.double(x %*% c(1, rep(0, p - 1L)))
    actual <- as.double(x %*% c(1, beta))
    set.seed(817)
    u <- runif(1, 0.7, 0.9)
    after <- .Random.seed
    expected <- u * EssDee(actual) * length(i)^(-1/5)
    set.seed(817)
    got <- .npindex_random_restart_bandwidth(x, beta, first, "fixed",
      length(i), controls, 0)
    expect_equal(got, expected, tolerance = 1e-14)
    expect_identical(.Random.seed, after)
    if (p > 1L)
      expect_false(isTRUE(all.equal(got, u * EssDee(first) * length(i)^(-1/5))))
    expect_error(.npindex_random_restart_bandwidth(x, beta, first, "fixed",
      length(i), controls, 1e6), "below the continuous scale-factor lower bound")
    for (type in c("generalized_nn", "adaptive_nn")) {
      set.seed(817)
      old <- .npindex_random_start_bandwidth(first, type, length(i), controls)
      after <- .Random.seed
      set.seed(817)
      expect_identical(.npindex_random_restart_bandwidth(x, beta, first, type,
        length(i), controls, 0), old)
      expect_identical(.Random.seed, after)
    }
  }
})

test_that("NOMAD preparation preserves first and non-bandwidth coordinates", {
  for (p in c(1L, 2L, 3L, 5L)) {
    i <- seq_len(48L)
    x <- vapply(seq_len(p), function(j) sin(i * sqrt(j + 1)) * j, numeric(48L))
    coord <- .npindex_beta_coordinate_setup(x)
    first.scale <- EssDee(x[,1]) * nrow(x)^(-1/5)
    starts <- matrix(0, 4L, p + 1L)
    if (p > 1L)
      starts[,seq_len(p-1L)] <- t(vapply(0:3, function(j)
        coord$to_search(rep(j/3, p-1L)), numeric(p-1L)))
    starts[,p] <- c(0.5, 0.7, 0.8, 0.9)
    starts[,p+1L] <- c(0,1,2,1)
    set.seed(817)
    before <- .Random.seed
    got <- .npindexbw_prepare_fixed_starts(starts,x,coord,first.scale)
    expect_identical(.Random.seed,before)
    expect_identical(got[1,],starts[1,])
    expect_identical(got[,-p,drop=FALSE],starts[,-p,drop=FALSE])
    for (j in 2:4) {
      beta <- if (p>1L) coord$to_public(starts[j,seq_len(p-1L)]) else numeric()
      scale <- EssDee(as.double(x %*% c(1,beta))) * nrow(x)^(-1/5)
      expect_equal(got[j,p] * first.scale, starts[j,p] * scale, tolerance=1e-14)
    }
    if (p==1L) expect_identical(got,starts)
  }
  x <- cbind(seq_len(8L),seq_len(8L))
  coord <- .npindex_beta_coordinate_setup(x)
  warnings <- character()
  withCallingHandlers(
    expect_error(.npindexbw_prepare_fixed_starts(rbind(c(0,.5,0),c(-1,.7,1)),
      x,coord,1), "nonpositive or nonfinite index scale"),
    warning=function(w) {warnings <<- c(warnings,conditionMessage(w)); invokeRestart("muffleWarning")})
  expect_length(warnings,2L)
  expect_match(warnings[1L],"variable 1 appears to be constant")
  expect_match(warnings[2L],"no non-missing arguments to min")
})

test_that("prepared NOMAD starts must respect the existing coordinate contract", {
  args <- list(engine="nomad",baseline_record=NULL,x0=c(.5,1),
    bbin=c(0L,1L),lb=c(0,0),ub=c(1,2),nmulti=2L,
    eval_fun=function(...) stop("unexpected objective evaluation"),
    build_payload=function(...) stop("unexpected payload construction"))
  for (prepare in list(function(s) s[,1,drop=FALSE],
                       function(s) {s[2,1]<-NA_real_;s},
                       function(s) {s[2,1]<- -1e-9;s},
                       function(s) {s[2,1]<-1+1e-9;s},
                       function(s) {s[2,2]<-.5;s})) {
    expect_error(do.call(.np_nomad_search,c(args,list(prepare_starts=prepare))),
      "prepared starts violate")
  }
})
