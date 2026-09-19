#ifndef NP_FANOUT_CONTROL_H
#define NP_FANOUT_CONTROL_H
#include <Rinternals.h>
#include <mpi.h>
/* Private channel lifecycle (policy lives in R/np.fanout.transaction.R).
 * base remains caller-owned. First begin participates collectively to create
 * a duplicated channel; native pool/request slots survive R unwinds. Master
 * passes an environment, preserved until successful finish; workers pass
 * R_NilValue. MPI requests use native slots, not borrowed R vector payloads.
 * finish requires completed control retirement and releases the R owner, but
 * retains the channel for reuse. close frees only an inactive, non-quarantined
 * pool; active/poisoned state is deliberately retained rather than freed under
 * pending MPI requests. close/finalize belong to the explicit MPI lifecycle,
 * not a GC finalizer. Native control-retirement checks are in fanout_control.c;
 * terminal result receipts remain the R transaction owner's responsibility. */
SEXP np_fanout_control_begin(MPI_Comm base, SEXP owner);
SEXP np_fanout_control_send(MPI_Comm base, int rank, int command);
SEXP np_fanout_control_poll(MPI_Comm base);
SEXP np_fanout_control_finish(MPI_Comm base);
SEXP np_fanout_control_owner(MPI_Comm base);
void np_fanout_control_close(MPI_Comm base);
void np_fanout_control_finalize(void);
#endif
