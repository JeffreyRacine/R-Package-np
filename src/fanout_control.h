#ifndef NP_FANOUT_CONTROL_H
#define NP_FANOUT_CONTROL_H
#include <Rinternals.h>
#include <mpi.h>
SEXP np_fanout_control_begin(MPI_Comm base, SEXP owner);
SEXP np_fanout_control_send(MPI_Comm base, int rank, int command);
SEXP np_fanout_control_poll(MPI_Comm base);
SEXP np_fanout_control_finish(MPI_Comm base);
SEXP np_fanout_control_owner(MPI_Comm base);
void np_fanout_control_close(MPI_Comm base);
void np_fanout_control_finalize(void);
#endif
