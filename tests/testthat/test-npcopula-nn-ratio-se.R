test_that("distributed NN copula errors retain ratio and request contracts", {
  if (!spawn_mpi_slaves()) skip("Could not initialize MPI context")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  ns <- asNamespace("npRmpi")
  owner <- get(".npcopula_density_se_nn", ns)
  x <- data.frame(x = sin(1:24) + (1:24)/24,
                  y = cos(1:24) + (1:24)/48)
  ex <- data.frame(x = c(.15, .65, .15), y = c(.2, .6, .2))
  n <- nrow(x); k <- 6L
  center <- function(z) z - mean(z)
  for (type in c("generalized_nn", "adaptive_nn")) {
    bw <- npudensbw(dat = x, bws = c(6, 6), bwtype = type,
                    bandwidth.compute = FALSE)
    ann <- identical(type, "adaptive_nn")
    radius <- pilot <- matrix(NA_real_, n, 2L)
    membership <- vector("list", 2L)
    if (ann) for (j in 1:2) {
      distance <- abs(outer(x[[j]], x[[j]], "-"))
      width <- min(floor((n - 1)^(2/3)), floor(min(k, n - 1 - k)/2))
      for (i in seq_len(n)) {
        d <- sort(distance[i, -i])
        radius[i, j] <- d[k]
        pilot[i, j] <- (2 * width/(n - 1))/(d[k + width] - d[k - width])
      }
      membership[[j]] <- distance <= radius[, j]
    }
    reference <- vapply(seq_len(nrow(ex)), function(q) {
      h <- if (ann) radius else matrix(rep(vapply(1:2, function(j)
        sort(abs(x[[j]] - ex[q, j]))[k], numeric(1L)), each = n), n, 2L)
      z <- sweep(as.matrix(x), 2L, as.numeric(ex[q, ]), function(a, b) b - a)/h
      w <- dnorm(z)/h
      M <- colMeans(w); D <- prod(M); K <- w[, 1L] * w[, 2L]
      J <- mean(K); phi <- K; marginal <- w
      if (ann) for (j in 1:2) {
        phi <- phi + as.numeric(crossprod(
          K * (1 - z[, j]^2)/(h[, j] * pilot[, j]), membership[[j]]))/n
        marginal[, j] <- marginal[, j] + as.numeric(crossprod(
          w[, j] * (1 - z[, j]^2)/(h[, j] * pilot[, j]), membership[[j]]))/n
      }
      score <- center(phi)/D - (J/D) *
        (center(marginal[, 1L])/M[1L] + center(marginal[, 2L])/M[2L])
      sqrt(sum(center(score)^2)) * if (ann) 1/sqrt(n * (n - 1)) else 1/n
    }, numeric(1L))
    answer <- owner(bw, x, ex)
    expect_equal(answer, reference, tolerance = 2e-10)
    expect_identical(answer[1L], answer[3L])
    expect_equal(owner(bw, x, ex[c(3, 2, 1), ]), answer[c(3, 2, 1)],
                  tolerance = 2e-10)
    single <- data.frame(x = ordered("a"))
    single.bw <- npudensbw(dat = single, bws = 0, bwtype = type,
                           bandwidth.compute = FALSE)
    expect_identical(owner(single.bw, single, single), 0)
    categorical <- data.frame(x = ordered(rep(letters[1:3], 4L)))
    categorical.bw <- npudensbw(dat = categorical, bws = .2, bwtype = type,
      ckertype = "uniform", bandwidth.compute = FALSE)
    expect_silent(categorical.se <- owner(categorical.bw, categorical, categorical))
    expect_identical(categorical.se, rep(0, nrow(categorical)))

    b1 <- npudensbw(dat = x[1L], bws = 6, bwtype = type,
                    bandwidth.compute = FALSE)
    expect_identical(owner(b1, x[1L], ex[1L]), rep(0, 3L))
    off <- npcopula(bws = bw, data = x)
    on <- npcopula(bws = bw, data = x, se = TRUE)
    expect_null(off$copulaerr)
    expect_identical(as.data.frame(off), as.data.frame(on))
    expect_identical(off$bws, on$bws)
    expect_identical(predict(on, se.fit = TRUE)$se.fit, se(on))
    set.seed(917)
    seed <- .Random.seed
    invisible(owner(bw, x, ex))
    expect_identical(.Random.seed, seed)
  }
})

test_that("NN ratio norm protects representable large and small errors", {
  norm <- get(".npcopula_density_se_nn_norm", asNamespace("npRmpi"))
  expect_equal(norm(c(-1e200, 1e200)), sqrt(2) * 1e200)
  expect_equal(norm(c(-1e-200, 1e-200)), sqrt(2) * 1e-200)
  expect_gt(norm(c(-1e-200, 1e-200)), 0)
  expect_identical(norm(c(4, 4)), 0)
  expect_identical(norm(c(NA_real_, 1)), NA_real_)
})

test_that("retained NN metadata preserves genuine kernel conditions", {
  ns <- asNamespace(getNamespaceName(environment(npcopula)))
  retained <- get(".npcopula_density_se_nn_retained", ns)
  expect_warning(retained(list(ckertype = "uniform"), warning("distinct condition")),
    "distinct condition")
  expect_error(retained(list(ckertype = "uniform"), stop("distinct failure")),
    "distinct failure")
  advisory <- unname(get(".np_io_prefix_text", ns)(
    "ignoring kernel order specified with uniform kernel type"))
  expect_warning(retained(list(ckertype = "gaussian"), warning(advisory)),
    "ignoring kernel order")
})


test_that("uniform ANN query factors reuse invocation-local training preparation", {
  if (exists("spawn_mpi_slaves", mode = "function")) {
    if (!spawn_mpi_slaves()) skip("Could not initialize MPI context")
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }
  ns <- asNamespace(getNamespaceName(environment(npcopula)))
  x <- data.frame(x = sin(1:32) + (1:32)/32,
                  y = cos(1:32) + (1:32)/64)
  bw <- npudensbw(dat = x, bws = c(6, 6), bwtype = "adaptive_nn",
                  ckertype = "uniform", bandwidth.compute = FALSE)
  retained <- get(".npcopula_density_se_nn_retained", ns)
  prepare <- get(".np_plot_density_ann_prepare", ns)
  factor <- get(".npcopula_density_se_nn_factor", ns)
  context <- retained(bw, prepare(bw, x, TRUE, uniform.faces = TRUE))
  before <- serialize(context, NULL)
  evaluation <- rbind(x[c(1, 5, 11), ], data.frame(x = 1e6, y = 1e6))
  for (i in seq_len(nrow(evaluation)))
    expect_identical(factor(bw, x, evaluation[i, , drop = FALSE], context),
                     factor(bw, x, evaluation[i, , drop = FALSE]))
  expect_identical(serialize(context, NULL), before)
  # Fresh preparation must continue to reject invalid training geometry.
  tied <- x; tied$x <- 1
  expect_error(prepare(bw, tied, TRUE), "nonregular ANN spacing geometry")
})


test_that("NN uncertainty work is absent on default and explicit OFF routes", {
  # Extended lane: actual namespace hooks (and MPI ranks) test work, not output.
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not initialize MPI context")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  ns <- asNamespace("npRmpi")
  manage <- function(action) {
    ns <- asNamespace("npRmpi")
    keys <- c(".npcopula_density_se_nn", ".npcopula_density_se_nn_local",
              ".npcopula_density_se_nn_factor", ".np_plot_density_ann_prepare",
              ".npRmpi_copula_density_se_nn")
    tag <- "np.nn.off.contract"
    if (action == "install") {
      originals <- stats::setNames(lapply(keys, function(k) get(k, ns)), keys)
      stopifnot(!any(vapply(originals, function(f) !is.null(attr(f, tag)), TRUE)))
      state <- new.env(parent = emptyenv())
      state$counts <- stats::setNames(integer(length(keys)), keys)
      state$originals <- originals
      for (key in keys) local({
        name <- key; original <- originals[[name]]
        hook <- function(...) {
          state$counts[name] <- state$counts[name] + 1L
          original(...)
        }
        attr(hook, tag) <- state
        unlockBinding(name, ns); assign(name, hook, ns); lockBinding(name, ns)
      })
    }
    state <- attr(get(keys[1L], ns), tag)
    stopifnot(is.environment(state))
    answer <- state$counts
    if (action == "restore") {
      for (key in keys) {
        unlockBinding(key, ns); assign(key, state$originals[[key]], ns)
        lockBinding(key, ns)
      }
      stopifnot(all(vapply(keys, function(k)
        identical(get(k, ns), state$originals[[k]]), TRUE)))
    }
    list(counts = answer)
  }
  # Do not transport the surrounding test environment to MPI workers.
  environment(manage) <- baseenv()
  remote <- get("mpi.remote.exec", ns)
  everywhere <- function(action) {
    call <- substitute(remote(F(A), simplify = FALSE), list(F = manage, A = action))
    lapply(c(list(manage(action)), eval(call)), `[[`, "counts")
  }
  # Uniform metadata re-enters the existing constructor advisory on public calls.
  # Match only that known condition; every other warning remains a test warning.
  advisory <- unname(get(".np_io_prefix_text", ns)(
    "ignoring kernel order specified with uniform kernel type"))
  quiet.uniform <- function(expr) withCallingHandlers(expr, warning = function(w) {
    if (identical(conditionMessage(w), advisory)) invokeRestart("muffleWarning")
  })
  installed <- FALSE
  on.exit(if (installed) everywhere("restore"), add = TRUE, after = FALSE)
  everywhere("install"); installed <- TRUE
  d <- data.frame(x = sin(1:32) + (1:32)/32,
                  y = cos(1:32) + (1:32)/64)
  grid <- data.frame(x = rep(c(.3, .5, .7), length.out = 65L),
                     y = rep(c(.35, .55, .65), length.out = 65L))
  for (type in c("generalized_nn", "adaptive_nn")) {
    bw <- npudensbw(dat = d, bws = c(6, 6), bwtype = type,
                    ckertype = "uniform", bandwidth.compute = FALSE)
    args <- list(bws = bw, data = d, u = grid, n.quasi.inv = 24)
    before <- everywhere("read")
    off <- quiet.uniform(do.call(npcopula, args))
    explicit <- quiet.uniform(do.call(npcopula, c(args, list(se = FALSE))))
    invisible(quiet.uniform(predict(off, u = grid, se.fit = FALSE, n.quasi.inv = 24)))
    expect_identical(everywhere("read"), before)
    expect_null(off$copulaerr)
    expect_null(explicit$copulaerr)
    # Positive controls prove the hooks are on the actual numerical path.
    invisible(quiet.uniform(do.call(npcopula, c(args, list(se = TRUE)))))
    after <- everywhere("read")
    delta <- Reduce(`+`, Map(`-`, after, before))
    expect_gt(unname(delta[".npcopula_density_se_nn"]), 0L)
    expect_gt(unname(delta[".npcopula_density_se_nn_local"]), 0L)
    expect_gt(unname(delta[".npcopula_density_se_nn_factor"]), 0L)
    if (type == "adaptive_nn") {
      expect_gt(unname(delta[".np_plot_density_ann_prepare"]), 0L)
      # 65 rows exceed the retained 64-row chunk: workers must do real work.
      if (!isTRUE(getOption("npRmpi.local.regression.mode", FALSE)))
        expect_gt(sum(vapply(Map(`-`, after[-1L], before[-1L]),
          function(x) unname(x[".np_plot_density_ann_prepare"]), 0)), 0)
    }
  }
  everywhere("restore"); installed <- FALSE
  expect_null(quiet.uniform(do.call(npcopula, args))$copulaerr)
})
