test_that("fixed KWS owner follows exact support geometry", {
  old_options <- options(np.messages = FALSE, np.largeh = FALSE,
                         np.macMseries.accelerate = FALSE)
  old_oracle <- Sys.getenv("NP_KWS_TREE_OWNER_ORACLE", unset = NA_character_)
  old_dimension_oracle <- Sys.getenv("NP_KWS_TREE_DIMENSION_ORACLE",
                                     unset = NA_character_)
  on.exit({
    options(old_options)
    if (is.na(old_oracle)) Sys.unsetenv("NP_KWS_TREE_OWNER_ORACLE")
    else Sys.setenv(NP_KWS_TREE_OWNER_ORACLE = old_oracle)
    if (is.na(old_dimension_oracle))
      Sys.unsetenv("NP_KWS_TREE_DIMENSION_ORACLE")
    else Sys.setenv(NP_KWS_TREE_DIMENSION_ORACLE = old_dimension_oracle)
  }, add = TRUE)

  set.seed(20260816L)
  n <- 71L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  ranges <- vapply(x, function(value) diff(range(value)), numeric(1L))
  full <- 1.25 * ranges / sqrt(5)
  prunable <- 0.20 * ranges / sqrt(5)

  owner <- function(..., ckertype = "epanechnikov") {
    options(np.tree = TRUE)
    Sys.setenv(NP_KWS_TREE_OWNER_ORACLE = "1",
               NP_KWS_TREE_DIMENSION_ORACLE = "1")
    messages <- capture.output(
      value <- npksum(txdat = x, ckertype = ckertype, ...),
      type = "message"
    )
    Sys.unsetenv("NP_KWS_TREE_OWNER_ORACLE")
    Sys.unsetenv("NP_KWS_TREE_DIMENSION_ORACLE")
    list(
      value = value,
      oracle = messages[grepl("NP_KWS_TREE_OWNER_ORACLE", messages,
                              fixed = TRUE)],
      dimension_oracle = messages[grepl("NP_KWS_TREE_DIMENSION_ORACLE",
                                        messages, fixed = TRUE)]
    )
  }

  full_normal <- owner(bws = full, leave.one.out = TRUE)
  expect_length(full_normal$oracle, 1L)
  expect_match(full_normal$oracle, "common=1 capability=1 owner=dense",
               fixed = TRUE)

  prunable_normal <- owner(bws = prunable, leave.one.out = TRUE)
  expect_length(prunable_normal$oracle, 1L)
  expect_match(prunable_normal$oracle, "common=1 capability=2 owner=tree",
               fixed = TRUE)
  expect_match(prunable_normal$dimension_oracle,
               "dimensions=2 inactive=0 active=2", fixed = TRUE)

  mixed_normal <- owner(bws = c(full[[1L]], prunable[[2L]]),
                        leave.one.out = TRUE)
  expect_match(mixed_normal$oracle,
               "common=1 capability=2 owner=tree", fixed = TRUE)
  expect_match(mixed_normal$dimension_oracle,
               "dimensions=2 inactive=1 active=1", fixed = TRUE)

  mixed_derivative <- owner(
    bws = c(full[[1L]], prunable[[2L]]),
    operator = "derivative",
    leave.one.out = TRUE
  )
  expect_length(mixed_derivative$dimension_oracle, 0L)

  unbounded_normal <- owner(bws = rep.int(0.05, 2L),
                            ckertype = "gaussian", leave.one.out = TRUE)
  expect_match(unbounded_normal$oracle,
               "common=1 capability=1 owner=dense", fixed = TRUE)

  full_integral <- owner(bws = full, operator = "integral",
                         leave.one.out = TRUE)
  expect_match(full_integral$oracle,
               "common=1 capability=1 owner=dense", fixed = TRUE)

  prunable_integral <- owner(bws = prunable, operator = "integral",
                             leave.one.out = TRUE)
  expect_match(prunable_integral$oracle,
               "common=1 capability=3 owner=tree", fixed = TRUE)

  convolution_only <- owner(
    bws = 0.30 * ranges,
    operator = "convolution",
    leave.one.out = TRUE
  )
  expect_match(convolution_only$oracle,
               "common=1 capability=1 owner=dense", fixed = TRUE)

  convolution_with_derivative <- owner(
    bws = 0.30 * ranges,
    operator = "convolution",
    permutation.operator = "derivative",
    leave.one.out = TRUE
  )
  expect_match(convolution_with_derivative$oracle,
               "common=1 capability=2 owner=tree", fixed = TRUE)

  separate_eval <- owner(bws = full, exdat = x[seq_len(n - 3L), ])
  expect_match(separate_eval$oracle,
               "requested=1 common=0 capability=0 owner=tree", fixed = TRUE)

  generalized <- owner(bws = rep.int(5, 2L), bwtype = "generalized_nn",
                       leave.one.out = TRUE)
  expect_match(generalized$oracle,
               "requested=1 common=0 capability=0 owner=tree", fixed = TRUE)

  adaptive <- owner(bws = rep.int(5, 2L), bwtype = "adaptive_nn",
                    leave.one.out = TRUE)
  expect_match(adaptive$oracle,
               "requested=1 common=0 capability=0 owner=tree", fixed = TRUE)
})

test_that("fixed KWS dense and tree siblings preserve all requested outputs", {
  old_options <- options(np.messages = FALSE, np.largeh = FALSE,
                         np.macMseries.accelerate = FALSE)
  on.exit(options(old_options), add = TRUE)

  set.seed(20260817L)
  n <- 61L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- cbind(1, sin(2 * pi * x$x1) + x$x2)
  weights <- cbind(1, seq_len(n) / n)
  ranges <- vapply(x, function(value) diff(range(value)), numeric(1L))
  full <- 1.25 * ranges / sqrt(5)
  prunable <- 0.20 * ranges / sqrt(5)
  cases <- list(
    normal_full = list(bws = full, operator = "normal"),
    normal_prunable = list(bws = prunable, operator = "normal"),
    normal_mixed = list(bws = c(full[[1L]], prunable[[2L]]),
                        operator = "normal"),
    convolution_full = list(bws = full, operator = "convolution"),
    derivative_full = list(bws = full, operator = "derivative"),
    integral_full = list(bws = full, operator = "integral"),
    union_prunable = list(
      bws = 0.30 * ranges,
      operator = "convolution",
      permutation.operator = "derivative",
      return.derivative.kernel.weights = TRUE
    )
  )

  for (case_name in names(cases)) {
    args <- c(
      list(
        txdat = x,
        tydat = y,
        weights = weights,
        ckertype = "epanechnikov",
        leave.one.out = TRUE,
        return.kernel.weights = TRUE
      ),
      cases[[case_name]]
    )
    options(np.tree = FALSE)
    dense <- do.call(npksum, args)
    options(np.tree = TRUE)
    tree <- do.call(npksum, args)

    for (field in c("ksum", "kw", "p.ksum", "p.kw")) {
      expect_equal(
        tree[[field]], dense[[field]], tolerance = 2e-11,
        info = paste(case_name, field)
      )
    }
  }
})

test_that("fixed tree ownership cannot depend on sample size", {
  header <- paste(
    readLines(np_test_source_path("src", "tree_capability.h"), warn = FALSE),
    collapse = "\n"
  )
  source <- paste(
    readLines(np_test_source_path("src", "jksum.c"), warn = FALSE),
    collapse = "\n"
  )
  selector <- sub(
    "(?s).*?(fixed_tree_common_geometry[[:space:]]*=[[:space:]]*np_ks_tree_use.*?const char \\* const owner_oracle).*",
    "\\1", source, perl = TRUE
  )

  expect_false(grepl(
    "\\b(num_obs|nobs|sample_size|sample_count)\\b",
    header, perl = TRUE
  ))
  expect_false(grepl(
    "num_obs_(train|eval)[[:space:]]*(<=|>=|<|>|==|!=)",
    selector, perl = TRUE
  ))
})

test_that("exact inactive tree dimensions do not trigger K(0) arithmetic", {
  source <- paste(
    readLines(np_test_source_path("src", "jksum.c"), warn = FALSE),
    collapse = "\n"
  )
  active_plan <- sub(
    "(?s).*?(const int tree_use_active_dims[[:space:]]*=.*?for \\(ii = 0; ii < p_nvar; ii\\+\\+\\)).*",
    "\\1", source, perl = TRUE
  )
  kernel_plan <- sub(
    "(?s).*?(const int use_largeh = any_cont_largeh.*?np_ckernelv\\().*",
    "\\1", source, perl = TRUE
  )

  expect_match(active_plan, "fixed_tree_nonpruning_dimension",
               fixed = TRUE)
  expect_match(active_plan, "cont_largeh_active", fixed = TRUE)
  expect_false(grepl("fixed_tree_nonpruning_dimension", kernel_plan,
                     fixed = TRUE))
})
