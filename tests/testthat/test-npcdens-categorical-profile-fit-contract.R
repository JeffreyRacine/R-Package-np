npc_categorical_profile_data <- function(n = 600L, seed = 1L,
                                         ordered_y = FALSE,
                                         ordered_x = FALSE) {
  set.seed(seed)
  y_raw <- sample(0:3, n, replace = TRUE)
  x1_raw <- sample(0:4, n, replace = TRUE)
  x2_raw <- sample(0:2, n, replace = TRUE)

  data.frame(
    y = if (ordered_y) ordered(y_raw) else factor(y_raw),
    x1 = if (ordered_x) ordered(x1_raw) else factor(x1_raw),
    x2 = factor(x2_raw)
  )
}

npc_categorical_profile_value <- function(fit) {
  if (!is.null(fit$condens)) fit$condens else fit$condist
}

npc_categorical_profile_se <- function(fit) {
  if (!is.null(fit$condens.se)) fit$condens.se else fit$condist.se
}

npc_categorical_profile_case <- function(fun, bwfun, dat, bws) {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE, np.categorical.compress = FALSE)
  on.exit(options(old_opts), add = TRUE)

  bw <- bwfun(
    xdat = dat[c("x1", "x2")],
    ydat = dat["y"],
    bws = bws,
    bandwidth.compute = FALSE
  )

  options(np.tree = FALSE, np.categorical.compress = FALSE)
  dense <- fun(bws = bw)

  options(np.tree = FALSE, np.categorical.compress = TRUE)
  profile <- fun(bws = bw)

  eval_rows <- c(1L, 7L, 19L, 101L, 211L)
  edat <- dat[eval_rows, , drop = FALSE]
  dense_pred <- predict(dense, exdat = edat[c("x1", "x2")], eydat = edat["y"])
  profile_pred <- predict(profile, exdat = edat[c("x1", "x2")], eydat = edat["y"])

  list(
    dense = dense,
    profile = profile,
    dense_pred = dense_pred,
    profile_pred = profile_pred
  )
}

test_that("all-categorical npcdens profile route preserves fitted values and prediction", {
  cases <- list(
    npc_categorical_profile_data(seed = 20260601L),
    npc_categorical_profile_data(seed = 20260602L, ordered_y = TRUE, ordered_x = TRUE),
    npc_categorical_profile_data(seed = 20260603L, ordered_x = TRUE)
  )

  for (dat in cases) {
    out <- npc_categorical_profile_case(npcdens, npcdensbw, dat, c(0.20, 0.25, 0.30))

    expect_equal(
      npc_categorical_profile_value(out$profile),
      npc_categorical_profile_value(out$dense),
      tolerance = 1e-12
    )
    expect_equal(
      npc_categorical_profile_se(out$profile),
      npc_categorical_profile_se(out$dense),
      tolerance = 1e-12
    )
    expect_equal(
      as.numeric(out$profile_pred),
      as.numeric(out$dense_pred),
      tolerance = 1e-12
    )
  }
})

test_that("all-categorical npcdist profile route preserves fitted values and prediction", {
  dat <- npc_categorical_profile_data(
    seed = 20260604L,
    ordered_y = TRUE,
    ordered_x = TRUE
  )

  out <- npc_categorical_profile_case(npcdist, npcdistbw, dat, c(0.20, 0.25, 0.30))

  expect_equal(
    npc_categorical_profile_value(out$profile),
    npc_categorical_profile_value(out$dense),
    tolerance = 1e-12
  )
  expect_equal(
    npc_categorical_profile_se(out$profile),
    npc_categorical_profile_se(out$dense),
    tolerance = 1e-12
  )
  expect_equal(
    as.numeric(out$profile_pred),
    as.numeric(out$dense_pred),
    tolerance = 1e-12
  )
})
