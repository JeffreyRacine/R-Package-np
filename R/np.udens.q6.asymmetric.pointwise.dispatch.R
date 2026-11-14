.np_udens_q6_asymmetric_pointwise_dispatch <- function(tasks,reference) {
  B <- sum(vapply(tasks,function(z)z$bsz,0L));out <- matrix(0,B,length(reference$grid)+1L)
  for(task in tasks) {
    ix <- seq.int(task$start,length.out=task$bsz)
    out[ix,] <- .np_udens_q6_asymmetric_pointwise_chunk(task,reference)
    .np_plot_progress_notify(stage="Constructing asymmetric Q6 marginal bands",done=max(ix),total=B,unit="draws")
  }
  out
}
