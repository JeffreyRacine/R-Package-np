locate_npcdens_degree_y_side_source <- function() {
  candidates <- c(
    file.path("src", "np.c"),
    file.path("..", "..", "src", "np.c")
  )
  hit <- candidates[file.exists(candidates)]
  if (!length(hit))
    return(NA_character_)
  normalizePath(hit[[1L]], mustWork = TRUE)
}

make_npcdens_degree_y_side_fixture <- function(n = 256L, seed = 42L) {
  set.seed(seed)
  x1 <- stats::runif(n, -1, 1)
  x2 <- stats::runif(n, -1, 1)
  signal <- 0.40 + 0.70 * x1 - 0.55 * x2 + 0.45 * x1^2 -
    0.30 * x2^2 + 0.65 * x1 * x2 + 0.55 * x1^3 - 0.45 * x2^3 +
    0.40 * x1^2 * x2 - 0.35 * x1 * x2^2
  # Preserve the campaign workload's shared RNG sequence before class draws.
  noise <- stats::rnorm(n, sd = 0.30)
  stopifnot(length(noise) == n)
  scores <- cbind(
    class0 = 0,
    class1 = 0.90 * signal,
    class2 = -0.75 * signal,
    class3 = 0.55 * signal - 0.35 * x1,
    class4 = -0.45 * signal + 0.40 * x2
  )
  shifted <- scores - apply(scores, 1L, max)
  probabilities <- exp(shifted)
  probabilities <- probabilities / rowSums(probabilities)
  uniforms <- stats::runif(n)
  response <- vapply(seq_len(n), function(i) {
    which(cumsum(probabilities[i, ]) >= uniforms[[i]])[[1L]] - 1L
  }, integer(1L))
  data.frame(
    y = factor(response, levels = 0:4, ordered = FALSE),
    x1 = x1,
    x2 = x2
  )
}

test_that("conditional-density degree search owns its future Y-side state", {
  path <- locate_npcdens_degree_y_side_source()
  skip_if(is.na(path), "source tree unavailable")
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")

  expect_match(
    text,
    "np_conditional_density_objective_needs_y_side\\(\\s*const int objective,\\s*const int degree_search\\)",
    perl = TRUE
  )
  expect_match(
    text,
    "objective == CBWM_CVML\\) &&\\s*\\(degree_search \\|\\|",
    perl = TRUE
  )
  expect_match(
    text,
    "np_conditional_density_objective_needs_y_side\\(\\s*ibwmfunc, degree_search\\)",
    perl = TRUE
  )
})

test_that("categorical CVML degree search can activate a positive-degree stream", {
  skip_on_cran()
  skip_if_not_installed("crs")
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  dat <- make_npcdens_degree_y_side_fixture()
  set.seed(42L)
  bw <- np::npcdensbw(
    y ~ x1 + x2,
    data = dat,
    nmulti = 1L,
    bwtype = "fixed",
    cxkertype = "epanechnikov",
    nomad = "auto",
    nomad.nmulti = 0L,
    random.seed = 42L
  )

  expect_identical(as.integer(bw$degree), c(2L, 1L))
  expect_gt(as.double(bw$fval), -394)
  expect_gt(as.integer(bw$num.feval.fast), 0L)
})
