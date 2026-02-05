/* Copyright (C) 2002 Hao Yu
 *  
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or   
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of 
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */
   
#include <R.h>
#include <R_ext/Rdynload.h>

#include "Rmpi.h"


int mpi_errhandler(int errcode)
{
	int errmsglen;
	char errmsg[MPI_MAX_ERROR_STRING];

/*	if (errcode != MPI_SUCCESS) {
		MPI_Error_string(errcode, errmsg, &errmsglen);
		error(errmsg);
	} */
	
	return errcode;
}

int erreturn(int errcode){
	if (errcode==MPI_SUCCESS)
		return 1;
	else
		return 0;
}


MPI_Datatype mpitype(SEXP sexp_type){
	MPI_Datatype datatype = MPI_DATATYPE_NULL;
	
	switch(INTEGER(sexp_type)[0]){
	case 1:
		datatype=MPI_INT;
		break;
	case 2:
		datatype=MPI_DOUBLE;
		break;
	case 3:
		datatype=MPI_CHAR;
		break;
	case 4:
                datatype=MPI_BYTE;
                break;
	}
	return datatype;
}

/*
void mystrcpy(char *new_str, char *old_str, int size) {
	int i;
	for (i=0; i< size; new_str[i]=old_str[i++]);
}
*/ 

/* Return a copy of a string using memory from R_alloc */

/*
char *acopy_string2(const char *in)
{
    char *out;
    int len = strlen(in);
    if (len > 0) {
        out = (char *) R_alloc(strlen(in), sizeof(char));
        strcpy(out, in);
    } else          
        out = "";  
    return out;
}

char *charsxp2char(SEXP x)
{
    return acopy_string2(CHAR(x));
}

*/

static const R_CallMethodDef callMethods[] = {
	/* In file "Rmpi.c". */
	{"mpi_initialize", (DL_FUNC) &mpi_initialize, 0},
	{"mpi_finalize", (DL_FUNC) &mpi_finalize, 0},
	{"mpi_get_processor_name", (DL_FUNC) &mpi_get_processor_name, 0},
	{"mpi_universe_size", (DL_FUNC) &mpi_universe_size, 0},
	{"mpi_any_source", (DL_FUNC) &mpi_any_source, 0},
	{"mpi_any_tag", (DL_FUNC) &mpi_any_tag, 0},
	{"mpi_undefined", (DL_FUNC) &mpi_undefined, 0},
	{"mpi_proc_null", (DL_FUNC) &mpi_proc_null, 0},
	{"mpi_info_create", (DL_FUNC) &mpi_info_create, 1},
	{"mpi_info_set", (DL_FUNC) &mpi_info_set, 3},
	{"mpi_info_get", (DL_FUNC) &mpi_info_get, 3},
	{"mpi_info_free", (DL_FUNC) &mpi_info_free, 1},
	{"mpi_realloc_comm", (DL_FUNC) &mpi_realloc_comm, 1},
	{"mpi_comm_maxsize", (DL_FUNC) &mpi_comm_maxsize, 0},
	{"mpi_realloc_status", (DL_FUNC) &mpi_realloc_status, 1},
	{"mpi_status_maxsize", (DL_FUNC) &mpi_status_maxsize, 0},
	{"mpi_realloc_request", (DL_FUNC) &mpi_realloc_request, 1},
	{"mpi_request_maxsize", (DL_FUNC) &mpi_request_maxsize, 0},
	{"mpi_realloc_datatype", (DL_FUNC) &mpi_realloc_datatype, 1},
	{"mpi_gather", (DL_FUNC) &mpi_gather, 5},
	{"mpi_gatherv", (DL_FUNC) &mpi_gatherv, 6},
	{"mpi_scatter", (DL_FUNC) &mpi_scatter, 5},
	{"mpi_scatterv", (DL_FUNC) &mpi_scatterv, 6},
	{"mpi_allgather", (DL_FUNC) &mpi_allgather, 4},
	{"mpi_allgatherv", (DL_FUNC) &mpi_allgatherv, 5},
	{"mpi_bcast", (DL_FUNC) &mpi_bcast, 5},
	{"mpi_send", (DL_FUNC) &mpi_send, 5},
	{"mpi_recv", (DL_FUNC) &mpi_recv, 6},
	{"mpi_reduce", (DL_FUNC) &mpi_reduce, 5},
	{"mpi_allreduce", (DL_FUNC) &mpi_allreduce, 4},
	{"mpi_iprobe", (DL_FUNC) &mpi_iprobe, 4},
	{"mpi_probe", (DL_FUNC) &mpi_probe, 4},
	{"mpi_get_count", (DL_FUNC) &mpi_get_count, 2},
	{"mpi_get_sourcetag", (DL_FUNC) &mpi_get_sourcetag, 1},
	{"mpi_barrier", (DL_FUNC) &mpi_barrier, 1},
	{"mpi_comm_is_null", (DL_FUNC) &mpi_comm_is_null, 1},
	{"mpi_comm_size", (DL_FUNC) &mpi_comm_size, 1},
	{"mpi_comm_rank", (DL_FUNC) &mpi_comm_rank, 1},
	{"mpi_comm_dup", (DL_FUNC) &mpi_comm_dup, 2},
	{"mpi_comm_c2f", (DL_FUNC) &mpi_comm_c2f, 1},
	{"mpi_comm_free", (DL_FUNC) &mpi_comm_free, 1},
	{"mpi_abort", (DL_FUNC) &mpi_abort, 1},
	{"mpi_comm_set_errhandler", (DL_FUNC) &mpi_comm_set_errhandler, 1},
	{"mpi_comm_test_inter", (DL_FUNC) &mpi_comm_test_inter, 1},
	{"mpi_comm_spawn", (DL_FUNC) &mpi_comm_spawn, 7},
	{"mpi_comm_get_parent", (DL_FUNC) &mpi_comm_get_parent, 1},
	{"mpi_is_master", (DL_FUNC) &mpi_is_master, 0},
	{"mpi_comm_disconnect", (DL_FUNC) &mpi_comm_disconnect, 1},
	{"mpi_intercomm_merge", (DL_FUNC) &mpi_intercomm_merge, 3},
	{"mpi_comm_remote_size", (DL_FUNC) &mpi_comm_remote_size, 1},
	{"mpi_sendrecv", (DL_FUNC) &mpi_sendrecv, 10},
	{"mpi_sendrecv_replace", (DL_FUNC) &mpi_sendrecv_replace, 8},
	{"mpi_cart_create", (DL_FUNC) &mpi_cart_create, 5},
	{"mpi_dims_create", (DL_FUNC) &mpi_dims_create, 3},
	{"mpi_cartdim_get", (DL_FUNC) &mpi_cartdim_get, 1},
	{"mpi_cart_get", (DL_FUNC) &mpi_cart_get, 2},
	{"mpi_cart_rank", (DL_FUNC) &mpi_cart_rank, 2},
	{"mpi_cart_coords", (DL_FUNC) &mpi_cart_coords, 3},
	{"mpi_cart_shift", (DL_FUNC) &mpi_cart_shift, 3},
	{"mpi_isend", (DL_FUNC) &mpi_isend, 6},
	{"mpi_irecv", (DL_FUNC) &mpi_irecv, 6},
	{"mpi_wait", (DL_FUNC) &mpi_wait, 2},
	{"mpi_test", (DL_FUNC) &mpi_test, 2},
	{"mpi_cancel", (DL_FUNC) &mpi_cancel, 1},
	{"mpi_test_cancelled", (DL_FUNC) &mpi_test_cancelled, 1},
	{"mpi_waitany", (DL_FUNC) &mpi_waitany, 2},
	{"mpi_testany", (DL_FUNC) &mpi_testany, 2},
	{"mpi_waitall", (DL_FUNC) &mpi_waitall, 1},
	{"mpi_testall", (DL_FUNC) &mpi_testall, 1},
	{"mpi_testsome", (DL_FUNC) &mpi_testsome, 1},
	{"mpi_waitsome", (DL_FUNC) &mpi_waitsome, 1},

	/* Finish R_CallMethodDef. */
	{NULL, NULL, 0}
};

void R_init_npRmpi(DllInfo *info){
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	R_useDynamicSymbols(info, TRUE);
}

/*
void *R_chk_calloc2(size_t nelem, size_t elsize)
{ 
    void *p;
#ifndef HAVE_WORKING_CALLOC
    if(nelem == 0)                  
        return(NULL);
#endif
    p = calloc(nelem, elsize);              
    if(!p) error(_("Calloc could not allocate (%d of %d) memory"),
                 nelem, elsize);
    return(p);
}

*/
