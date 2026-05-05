args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || length(args) > 2L) {
  stop("usage: Rscript summarize_cpu_monitor.R <cpu_monitor.csv> [out_csv]",
       call. = FALSE)
}

input <- args[[1L]]
output <- if (length(args) >= 2L) args[[2L]] else sub("[.]csv$", "_summary.csv", input)

d <- read.csv(input, stringsAsFactors = FALSE)
if (!nrow(d)) {
  write.csv(data.frame(), output, row.names = FALSE)
  quit(status = 0L)
}

d$pcpu <- suppressWarnings(as.numeric(d$pcpu))
d <- d[is.finite(d$pcpu) &
         (d$role == "worker" |
            grepl("^([^[:space:]]*/)?Rscript[[:space:]]|^/Library/Frameworks/R[.]framework/.*/Rscript[[:space:]]|^/Library/Frameworks/R[.]framework/.*/exec/R[[:space:]]|slavedaemon[.]R",
                  d$args)), , drop = FALSE]
if (!nrow(d)) {
  write.csv(data.frame(), output, row.names = FALSE)
  quit(status = 0L)
}

split_key <- paste(d$active_stage, d$timestamp, sep = "\r")
samples <- lapply(split(d, split_key), function(x) {
  data.frame(
    active_stage = x$active_stage[[1L]],
    timestamp = x$timestamp[[1L]],
    n_ranks = length(unique(x$pid)),
    min_pcpu = min(x$pcpu),
    median_pcpu = median(x$pcpu),
    max_pcpu = max(x$pcpu),
    stringsAsFactors = FALSE
  )
})
samples <- do.call(rbind, samples)
samples$active_sample <- samples$n_ranks >= 2L & samples$max_pcpu >= 50

by_stage <- lapply(split(samples, samples$active_stage), function(x) {
  active <- x[x$active_sample, , drop = FALSE]
  y <- if (nrow(active)) active else x
  data.frame(
    active_stage = x$active_stage[[1L]],
    samples = nrow(x),
    active_samples = nrow(active),
    max_ranks_seen = max(x$n_ranks),
    median_min_pcpu = median(y$min_pcpu),
    p10_min_pcpu = as.numeric(stats::quantile(y$min_pcpu, probs = 0.10,
                                               names = FALSE, type = 7)),
    median_rank_median_pcpu = median(y$median_pcpu),
    p10_rank_median_pcpu = as.numeric(stats::quantile(y$median_pcpu, probs = 0.10,
                                                       names = FALSE, type = 7)),
    suspicious_below_80 = nrow(active) > 0L && median(active$median_pcpu) < 80,
    suspicious_below_90 = nrow(active) > 0L && median(active$median_pcpu) < 90,
    stringsAsFactors = FALSE
  )
})
out <- do.call(rbind, by_stage)
out <- out[order(out$suspicious_below_80, out$suspicious_below_90,
                 out$median_rank_median_pcpu, decreasing = c(TRUE, TRUE, FALSE)), ]
write.csv(out, output, row.names = FALSE)
cat("Wrote", output, "\n")
