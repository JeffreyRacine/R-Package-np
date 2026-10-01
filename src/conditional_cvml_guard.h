#ifndef NP_CONDITIONAL_CVML_GUARD_H
#define NP_CONDITIONAL_CVML_GUARD_H

/* Publish an evaluation-level event, not a sum of replicated row counts.
 * The conditional CVML owner calls this after collective computation returns;
 * local-mode calls must never enter a communicator collective. */
static void np_conditional_cvml_guard_finish(const double before)
{
#ifdef MPI2
  if(iNum_Processors > 1 && !np_mpi_local_regression_active() &&
     comm != NULL && comm[1] != MPI_COMM_NULL){
    const int local_guarded = np_guarded_cvml_hits_get() > before;
    int any_guarded = 0;
    if(MPI_Allreduce(&local_guarded, &any_guarded, 1, MPI_INT, MPI_MAX,
                     comm[1]) != MPI_SUCCESS)
      error("conditional CVML guarded-evaluation synchronization failed");
    if(any_guarded && !local_guarded)
      np_guarded_cvml_hit();
  }
#else
  (void)before;
#endif
}

#endif
