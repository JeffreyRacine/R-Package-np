profile_empty_matrix <- function(n) {
  matrix(numeric(), nrow = n, ncol = 0L)
}

profile_unordered_kernel <- function(train, eval, lambda, ncat, kernel, operator) {
  same <- train == eval
  if (kernel == 0L) {
    normal <- ifelse(same, 1 - lambda, lambda / (ncat - 1))
    convolution <- ifelse(
      same,
      (1 - lambda)^2 + lambda^2 / (ncat - 1),
      lambda / (ncat - 1) *
        (2 * (1 - lambda) + (ncat - 2) * lambda / (ncat - 1))
    )
    derivative <- ifelse(same, -1, 1 / (ncat - 1))
  } else {
    inverse_norm <- 1 / ((ncat - 1) * lambda + 1)
    normal <- ifelse(same, 1, lambda) * inverse_norm
    convolution <- ifelse(
      same,
      1 + (ncat - 1) * lambda^2,
      lambda * (2 + (ncat - 2) * lambda)
    ) * inverse_norm^2
    derivative <- ifelse(
      same,
      (1 - ncat) * inverse_norm^2,
      inverse_norm * (lambda * (1 - ncat) * inverse_norm + 1)
    )
  }
  switch(
    as.character(operator),
    `0` = normal,
    `1` = convolution,
    `2` = derivative,
    `3` = normal,
    stop("unsupported unordered oracle operator")
  )
}

profile_ordered_normal <- function(train, eval, lambda, support, kernel) {
  distance <- abs(train - eval)
  if (kernel == 0L) {
    return(ifelse(distance == 0, 1 - lambda,
                  0.5 * (1 - lambda) * lambda^distance))
  }
  if (kernel == 1L) {
    return(lambda^distance)
  }
  if (kernel == 2L) {
    return(lambda^distance * (1 - lambda) / (1 + lambda))
  }
  denominator <- vapply(
    train,
    function(value) sum(lambda^abs(value - support)),
    numeric(1)
  )
  ifelse(denominator > 0, lambda^distance / denominator, 0)
}

profile_ordered_derivative <- function(train, eval, lambda, support, kernel) {
  distance <- as.integer(abs(train - eval))
  if (kernel == 0L) {
    return(ifelse(
      distance == 0L,
      -1,
      0.5 * lambda^(distance - 1L) * (distance - 2 * lambda)
    ))
  }
  if (kernel == 1L) {
    return(ifelse(distance == 0L, 0,
                  distance * lambda^(distance - 1L)))
  }
  if (kernel == 2L) {
    return(lambda^(distance - 1L) *
             (distance * (1 - lambda^2) - 2 * lambda))
  }

  vapply(seq_along(train), function(i) {
    distances <- as.integer(abs(train[[i]] - support))
    terms <- lambda^distances
    derivatives <- ifelse(
      distances == 0L,
      0,
      distances * lambda^(distances - 1L)
    )
    numerator <- lambda^distance[[i]]
    numerator_derivative <- if (distance[[i]] == 0L) {
      0
    } else {
      distance[[i]] * lambda^(distance[[i]] - 1L)
    }
    denominator <- sum(terms)
    (numerator_derivative * denominator -
       numerator * sum(derivatives)) / denominator^2
  }, numeric(1))
}

profile_ordered_convolution <- function(train, eval, lambda, support, kernel,
                                        finite_support) {
  if (finite_support) {
    return(vapply(seq_along(train), function(i) {
      sum(
        profile_ordered_normal(
          rep(train[[i]], length(support)),
          support,
          lambda,
          support,
          kernel
        ) *
          profile_ordered_normal(
            rep(eval, length(support)),
            support,
            lambda,
            support,
            kernel
          )
      )
    }, numeric(1)))
  }

  distance <- as.integer(abs(train - eval))
  if (kernel == 0L) {
    return(ifelse(
      distance == 0L,
      0.5 * (1 - lambda)^2 * (1 + 1 / (1 - lambda^2)),
      (0.5 * (1 - lambda))^2 * lambda^distance *
        (1 + distance + 2 / (1 - lambda^2))
    ))
  }
  if (kernel == 1L) {
    return(numeric(length(train)))
  }
  if (kernel == 2L) {
    scale <- ((1 - lambda) / (1 + lambda))^2
    return(scale * lambda^distance *
             ((1 + lambda^2) / (1 - lambda^2) + distance))
  }
  vapply(seq_along(train), function(i) {
    train_denominator <- sum(lambda^abs(train[[i]] - support))
    eval_denominator <- sum(lambda^abs(eval - support))
    sum(
      lambda^abs(train[[i]] - support) *
        lambda^abs(eval - support)
    ) / (train_denominator * eval_denominator)
  }, numeric(1))
}

profile_ordered_integral <- function(train, eval, lambda, support, kernel) {
  if (kernel == 0L) {
    distance <- as.integer(abs(eval - train))
    return(ifelse(
      eval == train,
      1 - 0.5 * lambda,
      ifelse(eval < train, 0.5 * lambda^distance, 1 - lambda^distance)
    ))
  }
  if (kernel == 1L) {
    return(vapply(train, function(value) {
      sum(lambda^abs(value - support[support <= eval]))
    }, numeric(1)))
  }
  if (kernel == 2L) {
    distance <- as.integer(abs(eval - train))
    geometric <- lambda^distance / (1 + lambda)
    return(ifelse(eval < train, geometric, 1 - lambda * geometric))
  }
  vapply(train, function(value) {
    denominator <- sum(lambda^abs(value - support))
    sum(lambda^abs(value - support[support <= eval])) / denominator
  }, numeric(1))
}

profile_ordered_kernel <- function(train, eval, lambda, support, kernel,
                                   operator, finite_support = FALSE) {
  switch(
    as.character(operator),
    `0` = profile_ordered_normal(train, eval, lambda, support, kernel),
    `1` = profile_ordered_convolution(
      train, eval, lambda, support, kernel, finite_support
    ),
    `2` = profile_ordered_derivative(train, eval, lambda, support, kernel),
    `3` = profile_ordered_integral(train, eval, lambda, support, kernel),
    stop("unsupported ordered oracle operator")
  )
}

profile_tile_call <- function(train_unordered, train_ordered,
                              eval_unordered, eval_ordered,
                              kernel_unordered, kernel_ordered,
                              operator, lambda, ncat, supports,
                              eval_start = 1L,
                              eval_count = nrow(eval_unordered)) {
  package <- if ("npRmpi" %in% loadedNamespaces()) "npRmpi" else "np"
  tile <- getFromNamespace(".np_categorical_profile_kernel_tile", package)
  tile(
    train_unordered,
    train_ordered,
    eval_unordered,
    eval_ordered,
    kernel_unordered,
    kernel_ordered,
    operator,
    lambda,
    ncat,
    supports,
    eval_start,
    eval_count
  )
}

test_that("categorical profile tiles match independent unordered formulas", {
  train <- matrix(c(1, 2, 1, 3, 2, 3), ncol = 1L)
  eval <- matrix(c(1, 3, 2), ncol = 1L)
  empty_train <- profile_empty_matrix(nrow(train))
  empty_eval <- profile_empty_matrix(nrow(eval))

  for (kernel in 0:1) {
    for (operator in 0:3) {
      actual <- profile_tile_call(
        train, empty_train, eval, empty_eval,
        kernel, integer(), operator, 0.27, 3L, list(NULL)
      )
      expected <- vapply(eval[, 1L], function(value) {
        profile_unordered_kernel(
          train[, 1L], value, 0.27, 3L, kernel, operator
        )
      }, numeric(nrow(train)))
      expect_equal(actual, expected, tolerance = 1e-14)
    }
  }
})

test_that("categorical profile tiles match independent ordered formulas", {
  support <- as.double(1:4)
  train <- matrix(c(1, 2, 4, 3, 1, 4), ncol = 1L)
  eval <- matrix(c(1, 3, 4), ncol = 1L)
  empty_train <- profile_empty_matrix(nrow(train))
  empty_eval <- profile_empty_matrix(nrow(eval))

  for (kernel in 0:3) {
    for (operator in 0:3) {
      actual <- profile_tile_call(
        empty_train, train, empty_eval, eval,
        integer(), kernel, operator, 0.31, 4L, list(support)
      )
      expected <- vapply(eval[, 1L], function(value) {
        profile_ordered_kernel(
          train[, 1L],
          value,
          0.31,
          support,
          kernel,
          operator,
          finite_support = operator == 1L && kernel == 1L
        )
      }, numeric(nrow(train)))
      expect_equal(actual, expected, tolerance = 2e-14)
    }
  }
})

test_that("mixed categorical dimensions retain product and tile order", {
  train_unordered <- cbind(
    c(1, 2, 1, 3, 2, 3),
    c(2, 2, 1, 1, 2, 1)
  )
  eval_unordered <- cbind(c(1, 3, 2, 1), c(2, 1, 2, 1))
  train_ordered <- cbind(
    c(1, 2, 4, 3, 1, 4),
    c(3, 2, 1, 2, 3, 1)
  )
  eval_ordered <- cbind(c(1, 3, 4, 2), c(2, 1, 3, 2))
  supports <- list(NULL, NULL, as.double(1:4), as.double(1:3))
  kernels_u <- c(0L, 1L)
  kernels_o <- c(1L, 3L)
  operators <- c(0L, 1L, 1L, 2L)
  lambdas <- c(0.21, 0.38, 0.29, 0.34)
  ncats <- c(3L, 2L, 4L, 3L)

  dense <- profile_tile_call(
    train_unordered, train_ordered,
    eval_unordered, eval_ordered,
    kernels_u, kernels_o, operators, lambdas, ncats, supports
  )
  tiled <- cbind(
    profile_tile_call(
      train_unordered, train_ordered,
      eval_unordered, eval_ordered,
      kernels_u, kernels_o, operators, lambdas, ncats, supports,
      eval_start = 1L, eval_count = 2L
    ),
    profile_tile_call(
      train_unordered, train_ordered,
      eval_unordered, eval_ordered,
      kernels_u, kernels_o, operators, lambdas, ncats, supports,
      eval_start = 3L, eval_count = 2L
    )
  )
  expected <- vapply(seq_len(nrow(eval_unordered)), function(row) {
    profile_unordered_kernel(
      train_unordered[, 1L], eval_unordered[row, 1L],
      lambdas[[1L]], ncats[[1L]], kernels_u[[1L]], operators[[1L]]
    ) *
      profile_unordered_kernel(
        train_unordered[, 2L], eval_unordered[row, 2L],
        lambdas[[2L]], ncats[[2L]], kernels_u[[2L]], operators[[2L]]
      ) *
      profile_ordered_kernel(
        train_ordered[, 1L], eval_ordered[row, 1L],
        lambdas[[3L]], supports[[3L]], kernels_o[[1L]], operators[[3L]],
        finite_support = TRUE
      ) *
      profile_ordered_kernel(
        train_ordered[, 2L], eval_ordered[row, 2L],
        lambdas[[4L]], supports[[4L]], kernels_o[[2L]], operators[[4L]],
        finite_support = TRUE
      )
  }, numeric(nrow(train_unordered)))

  expect_identical(tiled, dense)
  expect_equal(dense, expected, tolerance = 3e-14)
})

test_that("profile tiles are independent of categorical compression policy", {
  train_unordered <- matrix(c(1, 2, 1, 3, 2, 3), ncol = 1L)
  eval_unordered <- matrix(c(1, 3, 2), ncol = 1L)
  train_ordered <- matrix(c(1, 2, 4, 3, 1, 4), ncol = 1L)
  eval_ordered <- matrix(c(1, 3, 4), ncol = 1L)
  arguments <- list(
    train_unordered, train_ordered,
    eval_unordered, eval_ordered,
    1L, 3L, c(0L, 2L), c(0.27, 0.31),
    c(3L, 4L), list(NULL, as.double(1:4))
  )
  previous <- options(np.categorical.compress = TRUE)
  on.exit(options(previous), add = TRUE)
  compressed <- do.call(profile_tile_call, arguments)
  options(np.categorical.compress = FALSE)
  uncompressed <- do.call(profile_tile_call, arguments)

  expect_identical(compressed, uncompressed)
})

test_that("high-cardinality profiles stream through bounded tiles", {
  ntrain <- 20000L
  neval <- 257L
  tile_width <- 13L
  train_unordered <- matrix(as.double((seq_len(ntrain) - 1L) %% 17L + 1L),
                            ncol = 1L)
  eval_unordered <- matrix(as.double((seq_len(neval) * 7L) %% 17L + 1L),
                           ncol = 1L)
  train_ordered <- matrix(as.double((seq_len(ntrain) * 3L) %% 23L + 1L),
                          ncol = 1L)
  eval_ordered <- matrix(as.double((seq_len(neval) * 5L) %% 23L + 1L),
                         ncol = 1L)
  selected <- c(1L, ntrain %/% 2L, ntrain)
  support <- as.double(seq_len(23L))
  lambda <- c(0.19, 0.28)

  for (start in seq.int(1L, neval, by = tile_width)) {
    count <- min(tile_width, neval - start + 1L)
    actual <- profile_tile_call(
      train_unordered, train_ordered,
      eval_unordered, eval_ordered,
      0L, 2L, c(0L, 0L), lambda, c(17L, 23L),
      list(NULL, support), eval_start = start, eval_count = count
    )
    evaluation_rows <- seq.int(start, length.out = count)
    expected <- vapply(evaluation_rows, function(row) {
      profile_unordered_kernel(
        train_unordered[selected, 1L], eval_unordered[row, 1L],
        lambda[[1L]], 17L, 0L, 0L
      ) *
        profile_ordered_kernel(
          train_ordered[selected, 1L], eval_ordered[row, 1L],
          lambda[[2L]], support, 2L, 0L
        )
    }, numeric(length(selected)))

    expect_identical(dim(actual), c(ntrain, count))
    expect_lt(as.numeric(object.size(actual)), 3 * 1024^2)
    expect_equal(actual[selected, , drop = FALSE], expected, tolerance = 1e-14)
  }
})

test_that("categorical profile tile failures are explicit and bounded", {
  train <- matrix(as.double(rep(1:2, 5000L)), ncol = 1L)
  eval <- matrix(as.double(rep(1:2, 500L)), ncol = 1L)
  empty_train <- profile_empty_matrix(nrow(train))
  empty_eval <- profile_empty_matrix(nrow(eval))
  call <- function(...) {
    profile_tile_call(
      train, empty_train, eval, empty_eval,
      0L, integer(), 0L, 0.2, 2L, list(NULL), ...
    )
  }

  expect_error(call(eval_count = nrow(eval)), "64 MiB", fixed = TRUE)
  expect_error(call(eval_start = 0L, eval_count = 1L), "start must be positive")
  expect_error(call(eval_start = 1001L, eval_count = 1L), "range is invalid")
  expect_error(
    profile_tile_call(
      train, empty_train, eval, empty_eval,
      NA_integer_, integer(), 0L, 0.2, 2L, list(NULL),
      eval_count = 1L
    ),
    "unsupported categorical-profile kernel code"
  )
  expect_error(
    profile_tile_call(
      train, empty_train, eval, empty_eval,
      0L, integer(), 0L, NA_real_, 2L, list(NULL),
      eval_count = 1L
    ),
    "must be finite"
  )
  train[1L, 1L] <- NA_real_
  expect_error(call(eval_count = 1L), "must be finite")
})

test_that("categorical profile implementation is bounded and dormant", {
  candidates <- c(
    test_path("..", ".."),
    test_path("..", "..", ".."),
    Sys.getenv("R_PACKAGE_SOURCE", unset = "")
  )
  source_root <- candidates[
    file.exists(file.path(candidates, "src", "categorical_profile_tile.c"))
  ][1L]
  skip_if(is.na(source_root), "package source unavailable")

  native <- paste(
    readLines(
      file.path(source_root, "src", "categorical_profile_tile.c"),
      warn = FALSE
    ),
    collapse = "\n"
  )
  header <- paste(
    readLines(
      file.path(source_root, "src", "categorical_profile_tile.h"),
      warn = FALSE
    ),
    collapse = "\n"
  )
  production_r <- list.files(
    file.path(source_root, "R"),
    pattern = "\\.[Rr]$",
    full.names = TRUE
  )
  production_r <- production_r[
    basename(production_r) != "np.categorical.profile.tile.R"
  ]
  production_text <- paste(
    unlist(lapply(production_r, readLines, warn = FALSE)),
    collapse = "\n"
  )

  expect_match(header, "NP_CATEGORICAL_PROFILE_TILE_MAX_BYTES", fixed = TRUE)
  expect_match(header, "np_categorical_profile_spec_validate", fixed = TRUE)
  expect_match(
    header,
    "np_categorical_profile_tile_fill_prevalidated",
    fixed = TRUE
  )
  expect_match(native, "np_size_mul_checked", fixed = TRUE)
  expect_match(native, "np_size_array_bytes_checked", fixed = TRUE)
  expect_match(native, "double *output", fixed = TRUE)
  expect_false(grepl("malloc(", native, fixed = TRUE))
  expect_false(grepl(".np_categorical_profile_kernel_tile", production_text,
                     fixed = TRUE))

  jksum <- paste(
    readLines(file.path(source_root, "src", "jksum.c"), warn = FALSE),
    collapse = "\n"
  )
  expect_match(
    jksum,
    "kn_same = np_score_uaa(1, lambda, ncat);",
    fixed = TRUE
  )
  expect_match(
    jksum,
    "kn_same = np_score_unli_racine(1, lambda, ncat);",
    fixed = TRUE
  )
  expect_false(grepl(
    "const int kernel = (KERNEL >= 0 && KERNEL < nk) ? KERNEL : 0",
    jksum,
    fixed = TRUE
  ))
  expect_match(
    jksum,
    "Validate categorical operator dispatch once at the map-construction",
    fixed = TRUE
  )
  expect_match(
    jksum,
    "unsupported unordered kernel/operator combination",
    fixed = TRUE
  )
  expect_match(
    jksum,
    "unsupported ordered kernel/operator combination",
    fixed = TRUE
  )
})
