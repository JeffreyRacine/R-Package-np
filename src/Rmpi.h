#include <stdint.h>
#include <mpi.h>
#include <R.h>
#include <Rdefines.h>
#include <R_ext/Random.h>
#include <R_ext/Memory.h>
#include <Rinternals.h>

#define CHAR2(x) ((char *)CHAR(x))

/* #define COMM_MAXSIZE 10 */

int mpi_errhandler(int errcode);
int erreturn(int errcode);

MPI_Datatype mpitype(SEXP sexp_type); 

// void mystrcpy(char *new_str, char *old_str, int size);

SEXP AsInt (int n);

struct Dblint {
	double x;
	int rank;

};

// char *acopy_string2(const char *in);
// char *charsxp2char(SEXP x);

// #define CallocCharBuf2(n) 	(char *) R_chk_calloc2((size_t) ((n)+1), sizeof(char))
// extern void *R_chk_calloc2(size_t, size_t);

/* from Rmpi.c */
SEXP mpi_initialize(void);
SEXP mpi_finalize(void);
SEXP mpi_get_processor_name(void);
SEXP mpi_universe_size(void);
SEXP mpi_any_source(void);
SEXP mpi_any_tag(void);
SEXP mpi_undefined(void);
SEXP mpi_proc_null(void);
SEXP mpi_info_create(SEXP sexp_info);
SEXP mpi_info_set(SEXP sexp_info, SEXP sexp_key, SEXP sexp_value);
SEXP mpi_info_get(SEXP sexp_info, SEXP sexp_key, SEXP sexp_valuelen);
SEXP mpi_info_free(SEXP sexp_info);
SEXP mpi_realloc_comm(SEXP sexp_newncomm);
SEXP mpi_comm_maxsize(void);
SEXP mpi_realloc_status(SEXP sexp_newnstatus);
SEXP mpi_status_maxsize(void);
SEXP mpi_realloc_request(SEXP sexp_newnrequest);
SEXP mpi_request_maxsize(void);
SEXP mpi_realloc_datatype(SEXP sexp_newndatatype);
SEXP mpi_gather(SEXP sexp_sdata, SEXP sexp_type, SEXP sexp_rdata,
					SEXP sexp_root, SEXP sexp_comm);
SEXP mpi_gatherv(SEXP sexp_sdata, SEXP sexp_type, SEXP sexp_rdata, SEXP sexp_recvcounts,
					SEXP sexp_root, SEXP sexp_comm);
SEXP mpi_scatter(SEXP sexp_sdata, SEXP sexp_type, SEXP sexp_rdata,  SEXP sexp_root,
				   SEXP sexp_comm);
SEXP mpi_scatterv(SEXP sexp_sdata, SEXP sexp_sendcounts, SEXP sexp_type,
				  SEXP sexp_rdata, SEXP sexp_root, SEXP sexp_comm);
SEXP mpi_allgather(SEXP sexp_sdata, SEXP sexp_type, SEXP sexp_rdata, SEXP sexp_comm);
SEXP mpi_allgatherv(SEXP sexp_sdata, SEXP sexp_type, SEXP sexp_rdata, SEXP sexp_recvcounts,
				   SEXP sexp_comm);
SEXP mpi_bcast(SEXP sexp_data, SEXP sexp_type, SEXP	sexp_rank, SEXP sexp_comm,
			   SEXP sexp_buffunit);
SEXP mpi_send(SEXP sexp_data, SEXP sexp_type, SEXP sexp_dest, SEXP sexp_tag,
			  SEXP sexp_comm);
SEXP mpi_recv(SEXP sexp_data, SEXP sexp_type, SEXP sexp_source, SEXP sexp_tag,
			  SEXP sexp_comm, SEXP sexp_status);
SEXP mpi_reduce(SEXP sexp_send, SEXP sexp_type, SEXP sexp_op, SEXP sexp_dest, SEXP sexp_comm);
SEXP mpi_allreduce(SEXP sexp_send, SEXP sexp_type, SEXP sexp_op, SEXP sexp_comm);
SEXP mpi_iprobe(SEXP sexp_source, SEXP sexp_tag, SEXP sexp_comm, SEXP sexp_status);
SEXP mpi_probe(SEXP sexp_source, SEXP sexp_tag, SEXP sexp_comm, SEXP sexp_status);
SEXP mpi_get_count(SEXP sexp_status, SEXP sexp_type);
SEXP mpi_get_sourcetag (SEXP sexp_status);
SEXP mpi_barrier(SEXP sexp_comm);
SEXP mpi_comm_is_null(SEXP sexp_comm);
SEXP mpi_comm_size(SEXP sexp_comm);
SEXP mpi_comm_rank(SEXP sexp_comm);
SEXP mpi_comm_dup(SEXP sexp_comm, SEXP sexp_newcomm);
SEXP mpi_comm_c2f(SEXP sexp_comm);
SEXP mpi_comm_free(SEXP sexp_comm);
SEXP mpi_abort(SEXP sexp_comm);
SEXP mpi_comm_set_errhandler(SEXP sexp_comm);
SEXP mpi_comm_test_inter(SEXP sexp_comm);
SEXP mpi_comm_spawn (SEXP sexp_slave, SEXP sexp_argv, SEXP sexp_nslave, SEXP sexp_info,
					 SEXP sexp_root, SEXP sexp_intercomm, SEXP sexp_quiet);
SEXP mpi_comm_get_parent(SEXP sexp_comm);
SEXP mpi_is_master(void);
SEXP mpi_comm_disconnect(SEXP sexp_comm);
SEXP mpi_intercomm_merge(SEXP sexp_intercomm, SEXP sexp_high, SEXP sexp_comm);
SEXP mpi_comm_remote_size(SEXP sexp_comm);
SEXP mpi_sendrecv(SEXP sexp_senddata, SEXP sexp_sendtype, SEXP sexp_dest, SEXP sexp_sendtag,
        SEXP sexp_recvdata, SEXP sexp_recvtype, SEXP sexp_source, SEXP sexp_recvtag, SEXP sexp_comm,
        SEXP sexp_status);
SEXP mpi_sendrecv_replace(SEXP sexp_data, SEXP sexp_type, SEXP sexp_dest, SEXP sexp_sendtag,
        SEXP sexp_source, SEXP sexp_recvtag, SEXP sexp_comm, SEXP sexp_status);
SEXP mpi_cart_create(SEXP sexp_comm_old,  SEXP sexp_dims, SEXP sexp_periods, SEXP sexp_reorder,
           SEXP sexp_comm_cart);
SEXP mpi_dims_create(SEXP sexp_nnodes, SEXP sexp_ndims, SEXP sexp_dims);
SEXP mpi_cartdim_get(SEXP sexp_comm);
SEXP mpi_cart_get(SEXP sexp_comm, SEXP sexp_maxdims);
SEXP mpi_cart_rank(SEXP sexp_comm, SEXP sexp_coords);
SEXP mpi_cart_coords(SEXP sexp_comm, SEXP sexp_rank, SEXP sexp_maxdims);
SEXP mpi_cart_shift(SEXP sexp_comm, SEXP sexp_direction, SEXP sexp_disp);
SEXP mpi_isend(SEXP sexp_data, SEXP sexp_type, SEXP sexp_dest, SEXP sexp_tag,
                          SEXP sexp_comm, SEXP sexp_request);
SEXP mpi_irecv(SEXP sexp_data, SEXP sexp_type, SEXP sexp_source, SEXP sexp_tag,
                          SEXP sexp_comm, SEXP sexp_request);
SEXP mpi_wait(SEXP sexp_request, SEXP sexp_status);
SEXP mpi_test(SEXP sexp_request,  SEXP sexp_status);
SEXP mpi_cancel(SEXP sexp_request);
SEXP mpi_test_cancelled(SEXP sexp_status);
SEXP mpi_waitany(SEXP sexp_count, SEXP sexp_status);
SEXP mpi_testany(SEXP sexp_count, SEXP sexp_status);
SEXP mpi_waitall(SEXP sexp_count);
SEXP mpi_testall(SEXP sexp_count);
SEXP mpi_testsome(SEXP sexp_count);
SEXP mpi_waitsome(SEXP sexp_count);
