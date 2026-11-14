# New signed summary; primitive replay and incumbent transport remain unchanged.
.np_udens_q6_asymmetric_chunk <- function(task, reference) {
  Y <- .np_udens_q6_samples(task,reference$x,reference$h)
  z <- .np_udens_q6_asymmetric_maximum(Y,reference)
  cbind(plus=z$plus,minus=z$minus,refinement=rep(z$refinement,task$bsz))
}
.np_udens_q6_asymmetric_dispatch <- function(tasks, reference) {
  B <- sum(vapply(tasks,function(z)z$bsz,0L));out <- matrix(0,B,3L)
  for(task in tasks) {
    ix <- seq.int(task$start,length.out=task$bsz)
    out[ix,] <- .np_udens_q6_asymmetric_chunk(task,reference)
    .np_plot_progress_notify(stage="Constructing asymmetric Q6 simultaneous bands",
      done=max(ix),total=B,unit="draws")
  }
  out
}
