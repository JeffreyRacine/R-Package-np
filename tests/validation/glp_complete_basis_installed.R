args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L)
  stop("usage: glp_complete_basis_installed.R <library> <result.rds>")

.libPaths(c(args[1L], .libPaths()))
suppressPackageStartupMessages(library(np))

package <- "np"
nsfun <- function(name) get(name, envir = asNamespace(package), inherits = FALSE)
local_eval <- function(expr) eval(substitute(expr), envir = parent.frame())

build <- nsfun("npBuildLpTerms")
ncol_basis <- nsfun("npLpBasisNcol")
dim_bs <- nsfun("dimBS")
mypoly <- nsfun("mypoly")

checked <- 0L
for (p in 1:3) {
  degree_grid <- as.matrix(expand.grid(rep(list(0:4), p)))
  for (i in seq_len(nrow(degree_grid))) {
    degree <- as.integer(degree_grid[i, ])
    cart <- as.matrix(expand.grid(lapply(degree, function(d) 0:d)))
    expected <- cart[rowSums(cart) <= max(degree), , drop = FALSE]
    actual <- build(degree, basis = "glp")
    stopifnot(
      identical(unname(actual), unname(expected)),
      as.numeric(ncol_basis("glp", degree)) == nrow(expected),
      as.numeric(dim_bs("glp", kernel = TRUE, degree = degree,
                        segments = rep.int(1L, p))) == nrow(expected) - 1L
    )
    checked <- checked + 1L
  }
}

degree_141 <- c(1L, 4L, 1L)
cart_141 <- as.matrix(expand.grid(lapply(degree_141, function(d) 0:d)))
expected_141 <- cart_141[rowSums(cart_141) <= max(degree_141), , drop = FALSE]
stopifnot(identical(unname(build(degree_141, "glp")), unname(expected_141)))

train <- seq(-1.5, 2.5, length.out = 17L)
eval <- seq(-1.2, 2.2, length.out = 5L)
q <- 2 / diff(range(train))
d5 <- mypoly(x = train, ex = eval, degree = 5L,
             gradient.compute = TRUE, r = 5L,
             Bernstein = TRUE, complete.glp = TRUE)
expected_d5 <- sqrt(11) * (63 / 8) * factorial(5) * q^5
stopifnot(
  max(abs(d5[, 1:4, drop = FALSE])) == 0,
  max(abs(d5[, 5] - expected_d5)) < 2e-12
)

x_high <- seq(0.01, 0.99, length.out = 100L)
y_high <- x_high^13 / factorial(13)
bw_high <- local_eval(npregbw(
  xdat = data.frame(x = x_high),
  ydat = y_high,
  regtype = "lp",
  degree = 20L,
  degree.select = "manual",
  basis = "glp",
  bernstein.basis = TRUE,
  bws = 0.5,
  bandwidth.compute = FALSE
))
g_high <- local_eval(npreghat(
  bws = bw_high,
  txdat = data.frame(x = x_high),
  exdat = data.frame(x = seq(0.2, 0.8, length.out = 5L)),
  s = 13L
)) %*% y_high
high_error <- max(abs(as.numeric(g_high) - 1))
stopifnot(is.finite(high_error), high_error < 1e-3)

set.seed(20260717)
n <- 72L
tx <- data.frame(x1 = runif(n, -0.7, 0.8), x2 = runif(n, -0.6, 0.9))
y <- 0.2 + 0.4 * tx$x1 - 0.2 * tx$x1^2 + 0.1 * tx$x1^3 +
  0.3 * tx$x2 - 0.15 * tx$x2^2 + 0.5 * tx$x1 * tx$x2 +
  0.2 * tx$x1^2 * tx$x2 + 0.02 * cos(seq_len(n) * 1.3)
ex <- data.frame(x1 = seq(-0.5, 0.6, length.out = 7L),
                 x2 = seq(-0.4, 0.7, length.out = 7L))

fit_with <- function(bernstein) {
  bw <- local_eval(npregbw(
    xdat = tx, ydat = y, regtype = "lp", degree = c(3L, 2L),
    degree.select = "manual", basis = "glp",
    bernstein.basis = bernstein, bws = c(0.5, 0.5),
    bandwidth.compute = FALSE
  ))
  local_eval(npreg(
    txdat = tx, tydat = y, exdat = ex, bws = bw,
    gradients = TRUE, gradient.order = c(2L, 1L)
  ))
}

raw <- fit_with(FALSE)
shifted <- fit_with(TRUE)
errors <- c(
  mean = max(abs(raw$mean - shifted$mean)),
  grad = max(abs(raw$grad - shifted$grad)),
  merr = max(abs(raw$merr - shifted$merr)),
  gerr = max(abs(raw$gerr - shifted$gerr))
)
stopifnot(
  errors[["mean"]] < 1e-9,
  errors[["grad"]] < 1e-8,
  errors[["merr"]] < 1e-9,
  errors[["gerr"]] < 1e-8
)

result <- list(
  package = package,
  version = as.character(packageVersion(package)),
  vectors_checked = checked,
  degree_141_terms = nrow(expected_141),
  degree20_derivative13_error = high_error,
  degree32_functional_errors = errors
)
saveRDS(result, args[2L], version = 3L)
cat(sprintf("PASS: %s installed complete-GLP sentinel (%d degree vectors).\n",
            package, checked))
