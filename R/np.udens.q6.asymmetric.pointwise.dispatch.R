.np_udens_q6_asymmetric_pointwise_dispatch <- function(tasks,reference) {
  .npRmpi_bootstrap_fanout_enabled(n=length(reference$x),B=sum(vapply(tasks,function(z)z$bsz,0L)),
    chunk.size=32L,what="asymmetric Q6 density marginal process")
  worker <- function(task,reference).np_udens_q6_asymmetric_pointwise_chunk(task,reference)
  .npRmpi_bootstrap_run_fanout(tasks,worker,ncol.out=length(reference$grid)+1L,
    what="asymmetric Q6 density marginal process",profile.where="mpi.applyLB:q6-density-asymmetric-marginal",
    progress.label="Constructing asymmetric Q6 marginal bands",
    required.bindings=list(.np_udens_q6_asymmetric_pointwise_chunk=.np_udens_q6_asymmetric_pointwise_chunk,
      .np_udens_q6_samples=.np_udens_q6_samples,.np_udens_q6_density=.np_udens_q6_density,
      .np_udens_q6_curve=.np_udens_q6_curve,.np_udens_q6_convolution=.np_udens_q6_convolution,
      .np_udens_q6_rule=.np_udens_q6_rule,.np_udens_q6_kernel=.np_udens_q6_kernel,
      .np_udens_q6_check=.np_udens_q6_check),reference=reference)
}
