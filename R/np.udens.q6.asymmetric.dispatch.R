# New signed summary; primitive replay and incumbent transport remain unchanged.
.np_udens_q6_asymmetric_chunk <- function(task, reference) {
  Y <- .np_udens_q6_samples(task,reference$x,reference$h)
  z <- .np_udens_q6_asymmetric_maximum(Y,reference)
  cbind(plus=z$plus,minus=z$minus,refinement=rep(z$refinement,task$bsz))
}
.np_udens_q6_asymmetric_dispatch <- function(tasks, reference) {
  .npRmpi_bootstrap_fanout_enabled(n=length(reference$x),
    B=sum(vapply(tasks,function(z)z$bsz,0L)),chunk.size=32L,
    what="asymmetric Q6 density simultaneous process")
  worker <- function(task,reference) .np_udens_q6_asymmetric_chunk(task,reference)
  .npRmpi_bootstrap_run_fanout(tasks,worker,ncol.out=3L,
    what="asymmetric Q6 density simultaneous process",profile.where="mpi.applyLB:q6-density-asymmetric",
    progress.label="Constructing asymmetric Q6 simultaneous bands",
    required.bindings=list(.np_udens_q6_asymmetric_chunk=.np_udens_q6_asymmetric_chunk,
      .np_udens_q6_asymmetric_maximum=.np_udens_q6_asymmetric_maximum,
      .np_udens_q6_samples=.np_udens_q6_samples,.np_udens_q6_density=.np_udens_q6_density,
      .np_udens_q6_curve=.np_udens_q6_curve,.np_udens_q6_convolution=.np_udens_q6_convolution,
      .np_udens_q6_rule=.np_udens_q6_rule,.np_udens_q6_kernel=.np_udens_q6_kernel,
      .np_udens_q6_check=.np_udens_q6_check),reference=reference)
}
