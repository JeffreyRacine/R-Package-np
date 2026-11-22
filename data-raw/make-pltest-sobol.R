# Run at the package root in a fresh R session:
# Rscript --vanilla data-raw/make-pltest-sobol.R /path/to/LowDiscrepancy-sobol-orig1111.c
# Upstream data: randtoolbox 2.0.5, src/LowDiscrepancy-sobol-orig1111.c.
# Store 31 positive integer bits per generating column. Prefix i is obtained
# by XORing the columns selected by the Gray code of i, beginning at i = 1.
src <- paste(readLines(commandArgs(TRUE)[1L]), collapse = "\n")
array.values <- function(name) {
  pattern <- paste0("\\b", name, "\\[[^]]+\\] = \\{([^}]+)\\}")
  text <- regmatches(src, regexec(pattern, src))[[1L]][2L]
  as.integer(strsplit(gsub("[[:space:]]", "", text), ",", fixed = TRUE)[[1L]])
}
initial <- matrix(array.values("initmj"), ncol = 13L, byrow = TRUE)
degree <- array.values("mjshift")
a <- array.values("alla")
stopifnot(nrow(initial) == 1111L, length(degree) == 1111L, length(a) == 1111L)
words <- matrix(0L, 31L, 1111L)
words[, 1L] <- 2^(30:0)
storage.mode(words) <- "integer"
for (j in 2:1111) {
  s <- degree[j]
  for (bit in seq_len(s))
    words[bit, j] <- bitwShiftL(initial[j, bit], 31L - bit)
  for (bit in (s + 1L):31L) {
    prev <- words[bit - s, j]
    v <- bitwXor(prev, bitwShiftR(prev, s))
    for (k in seq_len(s - 1L))
      if (bitwAnd(bitwShiftR(a[j], s - 1L - k), 1L))
        v <- bitwXor(v, words[bit - k, j])
    words[bit, j] <- v
  }
}
# The fixed shift used in the mixed-data Monte Carlo. Storing it keeps
# projection selection independent of the caller's RNG kind and seed.
set.seed(514907L, kind = "Mersenne-Twister")
.pltest_sobol_table <- list(words = words, shift = runif(ncol(words)))
save(.pltest_sobol_table, file = "R/sysdata.rda", version = 2, compress = "xz")
