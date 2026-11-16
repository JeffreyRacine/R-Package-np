test_that('generated partially linear role transport preserves public terms ownership', {
 ns<-asNamespace('npRmpi');owner<-environment(get('npplreg.formula',ns))
 project<-get('.npRmpi_plreg_role_wire',owner);restore<-get('.npRmpi_plreg_role_restore',owner)
 env<-new.env(parent=baseenv());env$padding<-raw(800000)
 f<-~x;environment(f)<-env;d<-model.frame(f,data.frame(x=1:4))
 before<-serialize(d,NULL);p<-project(d);back<-restore(p,d)
 expect_identical(environment(attr(p,'terms')),emptyenv())
 expect_identical(environment(attr(d,'terms')),env)
 expect_identical(back,d);expect_identical(serialize(back,NULL),before)
 expect_identical(names(attributes(p)),names(attributes(d)))
 expect_identical(p$x,d$x);expect_identical(serialize(d,NULL),before)
 expect_lt(length(serialize(p,NULL)),length(before)-790000)
 expect_identical(project(data.frame(x=1:4)),data.frame(x=1:4))
 other<-p;attr(other,'terms')<-terms(~z)
 expect_identical(restore(other,d),other)
})


test_that("MPI formula evaluation keeps its original role terms", {
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves(1L), "MPI runtime not available")
  on.exit(close_mpi_slaves(), add = TRUE)
  owner <- new.env(parent = baseenv())
  i <- seq_len(24L)
  owner$d <- data.frame(y = cos(i * .3), x = sin(i * .7), z = cos(i * .4))
  form <- y ~ x | z
  environment(form) <- owner
  b <- npplregbw(form, data = owner$d, bws = matrix(.5, 2L, 1L),
                 bandwidth.compute = FALSE)
  fit <- npplreg(b, se = TRUE)
  expect_identical(environment(attr(fit$evalx, "terms", exact = TRUE)), owner)
  expect_identical(environment(fit$bws$formula), owner)
  newer <- npplreg(b, newdata = owner$d[1:5, ], se = FALSE)
  expect_identical(environment(attr(newer$evalx, "terms", exact = TRUE)), owner)
  external <- model.frame(~x, owner$d[1:5, ])
  attr(attr(external, "terms"), ".Environment") <- emptyenv()
  native <- npplreg(b, exdat = external, ezdat = owner$d[1:5, "z", drop = FALSE], se = FALSE)
  expect_identical(environment(attr(native$evalx, "terms", exact = TRUE)), emptyenv())
  expect_equal(fitted(native), fitted(newer), tolerance = 0)
})

test_that("formula constructor transfer does not grow with unrelated role bindings", {
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves(1L), "MPI runtime not available")
  on.exit(close_mpi_slaves(), add = TRUE)
  ns <- asNamespace("npRmpi")
  old <- options(l1.test.role.bytes = 0, l1.test.role.active = FALSE)
  on.exit(options(old), add = TRUE)
  for (name in c("mpi.send", "mpi.bcast"))
    trace(name, tracer = quote({
      if (isTRUE(getOption("l1.test.role.active")) && identical(as.integer(type), 4L))
        options(l1.test.role.bytes = getOption("l1.test.role.bytes") + length(x))
    }), where = ns, print = FALSE)
  on.exit({for (name in c("mpi.send", "mpi.bcast")) untrace(name, where = ns)}, add = TRUE)
  bytes <- numeric(2L)
  for (j in 1:2) {
    owner <- new.env(parent = baseenv())
    owner$padding <- raw(c(0L, 800000L)[j])
    i <- seq_len(24L)
    owner$d <- data.frame(y = cos(i*.3), x = sin(i*.7), z = cos(i*.4))
    form <- y ~ x | z
    environment(form) <- owner
    gc(); get(".npRmpi_lease_drain_pending", ns)()
    options(l1.test.role.bytes = 0, l1.test.role.active = TRUE)
    b <- npplregbw(form, data = owner$d, bws = matrix(.5, 2L, 1L),
                   bandwidth.compute = FALSE)
    options(l1.test.role.active = FALSE)
    bytes[j] <- getOption("l1.test.role.bytes")
    expect_identical(environment(b$formula), owner)
    rm(b)
  }
  # Allow bounded control-header variation; an 800KB environment must not travel.
  expect_lt(abs(diff(bytes)), 1024)
})

test_that("partially linear plot handoffs project both generated x roles", {
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves(1L), "MPI runtime not available")
  on.exit(close_mpi_slaves(), add = TRUE)
  ns <- asNamespace("npRmpi")
  owner <- new.env(parent = baseenv())
  i <- seq_len(24L)
  owner$d <- data.frame(y = cos(i*.3), x = sin(i*.7), z = cos(i*.4))
  form <- y ~ x | z
  environment(form) <- owner
  b <- npplregbw(form, data = owner$d, bws = matrix(.5, 2L, 1L),
                 bandwidth.compute = FALSE)
  dat <- get(".np_plreg_formula_training", ns)(b)
  ex <- dat$txdat[1:5, , drop = FALSE]
  ez <- dat$tzdat[1:5, , drop = FALSE]
  before <- serialize(list(dat = dat, ex = ex), NULL)
  bw.before <- get(".npRmpi_autodispatch_fingerprint", ns)(b)
  fit <- get(".np_plot_plreg_asymptotic_fit", ns)(
    b, dat$txdat, dat$tydat, dat$tzdat, ex, ez)
  expect_identical(environment(attr(fit$evalx, "terms", exact = TRUE)), owner)
  expect_identical(serialize(list(dat = dat, ex = ex), NULL), before)
  expect_identical(get(".npRmpi_autodispatch_fingerprint", ns)(b), bw.before)
  expect_identical(environment(b$formula), owner)
  old <- options(l1.test.plot.roles = list())
  on.exit(options(old), add = TRUE)
  trace(".npRmpi_bootstrap_prepare_worker", tracer = quote({
    if (all(c("txdat", "exdat") %in% names(required.bindings)))
      options(l1.test.plot.roles = c(getOption("l1.test.plot.roles"),
        list(vapply(required.bindings[c("txdat", "exdat")], function(x)
          identical(environment(attr(x, "terms", exact = TRUE)), emptyenv()),
          logical(1L)))))
  }), where = ns, print = FALSE)
  on.exit(untrace(".npRmpi_bootstrap_prepare_worker", where = ns), add = TRUE)
  counts <- matrix(1L, 24L, 9L)
  for (j in seq_len(9L)) { counts[j, j] <- 0L; counts[j + 9L, j] <- 2L }
  result <- get(".np_inid_boot_from_plreg_exact", ns)(dat$txdat, dat$tydat,
    dat$tzdat, ex, ez, b, B = 9L, counts = counts)
  roles <- getOption("l1.test.plot.roles")
  expect_gt(length(roles), 0L)
  expect_true(all(unlist(roles)))
  expect_identical(serialize(list(dat = dat, ex = ex), NULL), before)
  expect_identical(get(".npRmpi_autodispatch_fingerprint", ns)(b), bw.before)
  expect_identical(environment(b$formula), owner)
})
