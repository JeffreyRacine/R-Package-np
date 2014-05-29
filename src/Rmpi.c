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

#include "Rmpi.h"

#ifdef OPENMPI
#include <dlfcn.h>
#endif

static MPI_Comm	*comm;
static MPI_Status *status;
static MPI_Datatype *datatype;
static MPI_Info *info;
static MPI_Request *request;
static int COMM_MAXSIZE=10;
static int STATUS_MAXSIZE=2000;
static int REQUEST_MAXSIZE=2000;
static MPI_Datatype *xdouble;

#ifndef XLENGTH
#define XLENGTH LENGTH
#endif

SEXP mpidist(){
	int i=0;

#ifdef OPENMPI
	i=1;
#endif

#ifdef LAM
        i=2;
#endif

#ifdef MPICH
	i=3;
#endif

#if defined(MPICH2) || defined(INTELMPI)
	i=4;
#endif

	return AsInt(i);	
}

SEXP mpi_initialize(){
	int i,flag;
	MPI_Initialized(&flag);


#ifndef MPI2
        static int fake_argc = 1;
       	char *fake_argv[1];
        char *fake_argv0 = "R";
#endif

if (flag)
	return AsInt(1);
	else {

#ifndef __APPLE__
#ifdef OPENMPI
    if (!dlopen("libmpi.so.1", RTLD_GLOBAL | RTLD_LAZY) 
	&& !dlopen("libmpi.so.0", RTLD_GLOBAL | RTLD_LAZY)
	&& !dlopen("libmpi.so", RTLD_GLOBAL | RTLD_LAZY)) {
        Rprintf("%s\n",dlerror());
        return AsInt(0);
    }
#endif
#endif

#ifndef MPI2
   	fake_argv[0] = (char *)&fake_argv0;
       	MPI_Init(&fake_argc, (char ***)(void*)&fake_argv);
#else 
	MPI_Init((void *)0,(void *)0);
#endif 

		MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
		MPI_Comm_set_errhandler(MPI_COMM_SELF, MPI_ERRORS_RETURN);
		comm=(MPI_Comm *)Calloc(COMM_MAXSIZE, MPI_Comm); 
		status=(MPI_Status *)Calloc(STATUS_MAXSIZE, MPI_Status); 
		datatype=(MPI_Datatype *)Calloc(1, MPI_Datatype); 
		xdouble=(MPI_Datatype *)Calloc(1, MPI_Datatype); 
		info=(MPI_Info *)Calloc(1, MPI_Info);
		info[0]=MPI_INFO_NULL;
		request=(MPI_Request *)Calloc(REQUEST_MAXSIZE, MPI_Request);
		for (i=0; i< REQUEST_MAXSIZE; request[i++]=MPI_REQUEST_NULL);	
		comm[0]=MPI_COMM_WORLD;
		for (i=1;i < COMM_MAXSIZE; comm[i++]=MPI_COMM_NULL);

		return AsInt(1);
	} 
}

SEXP mpi_finalize(){
	MPI_Finalize();
	Free(comm);
	Free(status);
	Free(request);
	Free(datatype);
	Free(xdouble);
	Free(info);
	return AsInt(1);
}

SEXP mpi_get_processor_name (){
	int resultlen;
	char *name;
	SEXP sexp_name;
    PROTECT (sexp_name  = allocVector (STRSXP, 1));
	name = (char *)Calloc(MPI_MAX_PROCESSOR_NAME, char);
	MPI_Get_processor_name(name, &resultlen);
	SET_STRING_ELT(sexp_name, 0, mkChar(name));
	UNPROTECT(1);
	Free(name);

	return sexp_name;
}

/*
SEXP bin_nchar(SEXP sexp_data){
	return AsInt(LENGTH(STRING_ELT(sexp_data,0)));
}
*/

#ifdef MPI2
SEXP mpi_universe_size(){
	int *MPI_Universe_Size;
	int univ_flag;
	MPI_Comm_get_attr(comm[0], MPI_UNIVERSE_SIZE, &MPI_Universe_Size, &univ_flag);
	if (univ_flag)
    	return AsInt(*MPI_Universe_Size);
	else 
		return AsInt(0);
}
#endif

SEXP mpi_any_source(){
	return AsInt(MPI_ANY_SOURCE);
}

SEXP mpi_any_tag(){
	return AsInt(MPI_ANY_TAG);
}

SEXP mpi_undefined(){
	return AsInt(MPI_UNDEFINED);
}

SEXP mpi_proc_null(){
	return AsInt(MPI_PROC_NULL);
}

SEXP mpi_info_create(SEXP sexp_info){
	return AsInt(erreturn(mpi_errhandler(MPI_Info_create( &info[INTEGER(sexp_info)[0]]))));
}

SEXP mpi_info_set(SEXP sexp_info, SEXP sexp_key, SEXP sexp_value){
	return AsInt(erreturn(mpi_errhandler(MPI_Info_set(info[INTEGER(sexp_info)[0]],
		CHAR2( STRING_ELT (sexp_key,0)), CHAR2(STRING_ELT(sexp_value,0))))));
}

SEXP mpi_info_get(SEXP sexp_info, SEXP sexp_key, SEXP sexp_valuelen){
	int flag;
        char *value;
	SEXP sexp_value;

    	PROTECT (sexp_value  = allocVector (STRSXP, 1));
 	value = (char *)Calloc(INTEGER(sexp_valuelen)[0], char);  
	mpi_errhandler(MPI_Info_get(info[INTEGER(sexp_info)[0]], 
		CHAR2( STRING_ELT (sexp_key,0)), 
		INTEGER(sexp_valuelen)[0], value, &flag));
        SET_STRING_ELT(sexp_value, 0, mkChar(value));
        UNPROTECT(1); 
	Free(value);
	return sexp_value;
}

SEXP mpi_info_free(SEXP sexp_info){
	return AsInt(erreturn(mpi_errhandler(MPI_Info_free( &info[INTEGER(sexp_info)[0]]))));
}

SEXP mpi_realloc_comm(SEXP sexp_newncomm){
	int i, newcomm=INTEGER(sexp_newncomm)[0];
	if (newcomm > COMM_MAXSIZE){
		comm=(MPI_Comm *)Realloc(comm, newcomm, MPI_Comm); 	
		for (i=COMM_MAXSIZE; i < newcomm; comm[i++]=MPI_COMM_NULL);
		COMM_MAXSIZE=newcomm;
	}
	return AsInt(1);
}

SEXP mpi_comm_maxsize(){
	return AsInt(COMM_MAXSIZE);
}

SEXP mpi_realloc_status(SEXP sexp_newnstatus){
	int newsize=INTEGER(sexp_newnstatus)[0];
	if (newsize > STATUS_MAXSIZE){
		status=(MPI_Status *)Realloc(status, newsize, MPI_Status); 
		STATUS_MAXSIZE=newsize;
	}
	return AsInt(1);
}

SEXP mpi_status_maxsize(){
	return AsInt(STATUS_MAXSIZE);
}

SEXP mpi_realloc_request(SEXP sexp_newnrequest){
	int i, newsize=INTEGER(sexp_newnrequest)[0];
	if (newsize > REQUEST_MAXSIZE){
		request=(MPI_Request *)Realloc(request, newsize , MPI_Request); 
		for (i=REQUEST_MAXSIZE; i< newsize; request[i++]=MPI_REQUEST_NULL);	
		REQUEST_MAXSIZE=newsize;
	}
	return AsInt(1);
}

SEXP mpi_request_maxsize(){
	return AsInt(REQUEST_MAXSIZE);
}

SEXP mpi_realloc_datatype(SEXP sexp_newndatatype){
	datatype=(MPI_Datatype *)Realloc(datatype, INTEGER(sexp_newndatatype)[0], MPI_Datatype); 
	return AsInt(1);
}

/******************** Collective ***************************************/
SEXP mpi_gather(SEXP sexp_sdata,
				   SEXP sexp_type,
				   SEXP sexp_rdata,
				   SEXP sexp_root,
				   SEXP sexp_comm){
	int len, rlen, commn=INTEGER(sexp_comm)[0], root=INTEGER(sexp_root)[0];
	char *rdata;
	SEXP sexp_rdata2 = NULL;

	switch (INTEGER(sexp_type)[0]){
	case 1:
 		len=LENGTH(sexp_sdata);
		mpi_errhandler(MPI_Gather(INTEGER(sexp_sdata), len, MPI_INT, 
			INTEGER(sexp_rdata), len, MPI_INT, root, comm[commn]));
		break;
	case 2:
 		len=LENGTH(sexp_sdata);
		mpi_errhandler(MPI_Gather(REAL(sexp_sdata), len, MPI_DOUBLE, 
			REAL(sexp_rdata), len, MPI_DOUBLE, root, comm[commn]));
		break;
	case 3: 
		len=LENGTH(STRING_ELT(sexp_sdata,0));
		rlen=LENGTH(STRING_ELT(sexp_rdata,0));

        	PROTECT (sexp_rdata2  = allocVector (STRSXP, 1));
        	rdata = (char *)Calloc(rlen, char);
                MPI_Gather(CHAR2 (STRING_ELT ((sexp_sdata),0)), len, MPI_CHAR,
                        rdata, len, MPI_CHAR, root, comm[commn]);
        	SET_STRING_ELT(sexp_rdata2, 0, mkChar(rdata));
        	UNPROTECT(1);
		Free(rdata);
		break;
	case 4:
                len=LENGTH(sexp_sdata);
                mpi_errhandler(MPI_Gather(RAW(sexp_sdata), len, MPI_BYTE,
                        RAW(sexp_rdata), len, MPI_BYTE, root, comm[commn]));
                break;
	default:
		PROTECT(sexp_sdata=AS_NUMERIC(sexp_sdata));
		mpi_errhandler(MPI_Bcast(REAL(sexp_sdata), 1, datatype[0], root, comm[commn]));
		UNPROTECT(1);
		break;		
	}

	if (INTEGER(sexp_type)[0]==3)
		return sexp_rdata2;
	else
		return sexp_rdata;
}

SEXP mpi_gatherv(SEXP sexp_sdata,
				   SEXP sexp_type,
				   SEXP sexp_rdata,
				   SEXP sexp_recvcounts,
				   SEXP sexp_root,
				   SEXP sexp_comm){
	int len, rlen, commn=INTEGER(sexp_comm)[0], root=INTEGER(sexp_root)[0];
	int *displs=NULL, gsize, rank, i;
	char *rdata;
	SEXP sexp_rdata2 = NULL;
	
	MPI_Comm_size(comm[commn], &gsize);
	MPI_Comm_rank(comm[commn], &rank);
	if (rank==root){
		displs=(int *)Calloc(gsize, int);
		displs[0]=0;
		for (i=1; i < gsize; i++)
			displs[i]=displs[i-1]+INTEGER(sexp_recvcounts)[i-1];
	}

	switch (INTEGER(sexp_type)[0]){
	case 1:
		len=LENGTH(sexp_sdata);
		mpi_errhandler(MPI_Gatherv(INTEGER(sexp_sdata), len, MPI_INT, 
			INTEGER(sexp_rdata), INTEGER(sexp_recvcounts), displs, MPI_INT, 
				root, comm[commn]));
		break;
	case 2:
		len=LENGTH(sexp_sdata);
		mpi_errhandler(MPI_Gatherv(REAL(sexp_sdata), len, MPI_DOUBLE, 
			REAL(sexp_rdata), INTEGER(sexp_recvcounts), displs, 
				MPI_DOUBLE, root, comm[commn]));
		break;
	case 3:
		len=LENGTH(STRING_ELT(sexp_sdata,0));	
		rlen=LENGTH(STRING_ELT(sexp_rdata,0));
        
                PROTECT (sexp_rdata2  = allocVector (STRSXP, 1));
                rdata = (char *)Calloc(rlen, char);
                MPI_Gatherv(CHAR2 (STRING_ELT ((sexp_sdata),0)),len,MPI_CHAR,
		     rdata, INTEGER(sexp_recvcounts),
		     displs, MPI_CHAR, root, comm[commn]);
                SET_STRING_ELT(sexp_rdata2, 0, mkChar(rdata));
                UNPROTECT(1);
		Free(rdata);
		break;
	case 4:
                len=LENGTH(sexp_sdata);
                mpi_errhandler(MPI_Gatherv(RAW(sexp_sdata), len, MPI_BYTE,
                        RAW(sexp_rdata), INTEGER(sexp_recvcounts), displs,
                                MPI_BYTE, root, comm[commn]));
                break;

	default:
		PROTECT(sexp_sdata=AS_NUMERIC(sexp_sdata));
		mpi_errhandler(MPI_Bcast(REAL(sexp_sdata), 1, datatype[0], rank, comm[commn]));
		UNPROTECT(1);
		break;	
	}
	if (rank == root)
		Free(displs);

	if (INTEGER(sexp_type)[0]==3)
		return sexp_rdata2;
	else
		return sexp_rdata;
}

SEXP mpi_scatter(SEXP sexp_sdata,
				   SEXP sexp_type,
				   SEXP sexp_rdata,
				   SEXP sexp_root,
				   SEXP sexp_comm){
	int 	len, rlen;
	int	commn=INTEGER(sexp_comm)[0], root=INTEGER(sexp_root)[0];
	char 	*rdata;
	SEXP 	sexp_rdata2 = NULL;

	switch (INTEGER(sexp_type)[0]){
	case 1:
 		len=LENGTH(sexp_rdata);
		mpi_errhandler(MPI_Scatter(INTEGER(sexp_sdata), len, MPI_INT, 
			INTEGER(sexp_rdata), len, MPI_INT, root, comm[commn]));
		break;
	case 2:
		len=LENGTH(sexp_rdata);
		mpi_errhandler(MPI_Scatter(REAL(sexp_sdata), len, MPI_DOUBLE, 
			REAL(sexp_rdata), len, MPI_DOUBLE, root, comm[commn]));
		break;
	case 3:
 		len=LENGTH(STRING_ELT(sexp_rdata,0));
                rlen=LENGTH(STRING_ELT(sexp_rdata,0));

                PROTECT (sexp_rdata2  = allocVector (STRSXP, 1));
                rdata = (char *)Calloc(rlen, char);
		MPI_Scatter(CHAR2(STRING_ELT ((sexp_sdata),0)), len, MPI_CHAR,
                       rdata, len, MPI_CHAR, root, comm[commn]);
                SET_STRING_ELT(sexp_rdata2, 0, mkChar(rdata));
                UNPROTECT(1);
		Free(rdata);
		break;
	case 4:
                len=LENGTH(sexp_rdata);
                mpi_errhandler(MPI_Scatter(RAW(sexp_sdata), len, MPI_BYTE,
                        RAW(sexp_rdata), len, MPI_BYTE, root, comm[commn]));
                break;

	default:
		PROTECT(sexp_sdata=AS_NUMERIC(sexp_sdata));
		mpi_errhandler(MPI_Bcast(REAL(sexp_sdata), 1, datatype[0], root, comm[commn]));
		UNPROTECT(1);
		break;		
	}
        if (INTEGER(sexp_type)[0]==3)
                return sexp_rdata2;
        else
                return sexp_rdata;
}

SEXP mpi_scatterv(SEXP sexp_sdata,
				  SEXP sexp_sendcounts,
				  SEXP sexp_type,
				  SEXP sexp_rdata,
				  SEXP sexp_root,
				  SEXP sexp_comm){
	int len, rlen, commn=INTEGER(sexp_comm)[0], root=INTEGER(sexp_root)[0];
	int gsize,rank,i,*displs=NULL;
    	char *rdata;
	SEXP sexp_rdata2 = NULL;

	MPI_Comm_size(comm[commn], &gsize);
	MPI_Comm_rank(comm[commn], &rank);
	if (rank==root){
		displs=(int *)Calloc(gsize, int);
		displs[0]=0;
		for (i=1; i < gsize; i++)
			displs[i]=displs[i-1]+INTEGER(sexp_sendcounts)[i-1];
	}
	
	switch (INTEGER(sexp_type)[0]){
	case 1:
		 len=LENGTH(sexp_rdata);
		mpi_errhandler(MPI_Scatterv(INTEGER(sexp_sdata), INTEGER(sexp_sendcounts),
			displs, MPI_INT, INTEGER(sexp_rdata), len, MPI_INT, 
				root, comm[commn]));
		break;
	case 2:
		 len=LENGTH(sexp_rdata);
		mpi_errhandler(MPI_Scatterv(REAL(sexp_sdata), INTEGER(sexp_sendcounts),
			displs, MPI_DOUBLE, REAL(sexp_rdata), len,  
				MPI_DOUBLE, root, comm[commn]));
		break;
	case 3:
                len=LENGTH(STRING_ELT(sexp_rdata,0));
                rlen=LENGTH(STRING_ELT(sexp_rdata,0));

                PROTECT (sexp_rdata2  = allocVector (STRSXP, 1));
                // rdata = (char *)R_alloc(rlen, sizeof(char));
                rdata = (char *)Calloc(rlen, char);
                MPI_Scatterv(CHAR2 (STRING_ELT ((sexp_sdata),0)), INTEGER(sexp_sendcounts),displs, 
			MPI_CHAR,rdata, len, MPI_CHAR, root, comm[commn]);
                SET_STRING_ELT(sexp_rdata2, 0, mkChar(rdata));
                UNPROTECT(1);
		Free(rdata);
		break;
	case 4:
                len=LENGTH(sexp_rdata);
                mpi_errhandler(MPI_Scatterv(RAW(sexp_sdata), INTEGER(sexp_sendcounts),
                        displs, MPI_BYTE, RAW(sexp_rdata), len,
                                MPI_BYTE, root, comm[commn]));
                break;

	default:
		PROTECT(sexp_sdata=AS_NUMERIC(sexp_sdata));
		mpi_errhandler(MPI_Bcast(REAL(sexp_sdata), 1, datatype[0], rank, comm[commn]));
		UNPROTECT(1);
		break;		
	}
	if (rank == root)

		Free(displs);

        if (INTEGER(sexp_type)[0]==3)
                return sexp_rdata2;
        else
                return sexp_rdata;
}

SEXP mpi_allgather(SEXP sexp_sdata,
				   SEXP sexp_type,
				   SEXP sexp_rdata,
				   SEXP sexp_comm){
	int len, rlen, commn=INTEGER(sexp_comm)[0];
	char *rdata;
	SEXP sexp_rdata2 = NULL;
	
	switch (INTEGER(sexp_type)[0]){
	case 1:
		len=LENGTH(sexp_sdata);
		mpi_errhandler(MPI_Allgather(INTEGER(sexp_sdata), len, MPI_INT, 
			INTEGER(sexp_rdata), len, MPI_INT, comm[commn]));
		break;
	case 2:
		len=LENGTH(sexp_sdata);
		mpi_errhandler(MPI_Allgather(REAL(sexp_sdata), len, MPI_DOUBLE, 
			REAL(sexp_rdata), len, MPI_DOUBLE, comm[commn]));
		break;
	case 3:
	 	len=LENGTH(STRING_ELT(sexp_sdata,0));
                rlen=LENGTH(STRING_ELT(sexp_rdata,0));

                PROTECT (sexp_rdata2  = allocVector (STRSXP, 1));
                rdata = (char *)Calloc(rlen, char);
                MPI_Allgather(CHAR2 (STRING_ELT ((sexp_sdata),0)),len,
			MPI_CHAR,rdata, len, MPI_CHAR, comm[commn]);
                SET_STRING_ELT(sexp_rdata2, 0, mkChar(rdata));
                UNPROTECT(1);
		Free(rdata);
		break;
 	case 4:
                len=LENGTH(sexp_sdata);
                mpi_errhandler(MPI_Allgather(RAW(sexp_sdata), len, MPI_BYTE,
                        RAW(sexp_rdata), len, MPI_BYTE, comm[commn]));
                break;

	default:
		PROTECT(sexp_sdata=AS_NUMERIC(sexp_sdata));
		mpi_errhandler(MPI_Bcast(REAL(sexp_sdata), 1, datatype[0], 0, comm[commn]));
		UNPROTECT(1);
		break;
	}

        if (INTEGER(sexp_type)[0]==3)
                return sexp_rdata2;
        else
                return sexp_rdata;
}

SEXP mpi_allgatherv(SEXP sexp_sdata,
				   SEXP sexp_type,
				   SEXP sexp_rdata,
				   SEXP sexp_recvcounts,
				   SEXP sexp_comm){
	int len, rlen, commn=INTEGER(sexp_comm)[0], *displs, gsize, i;
	char *rdata;
	SEXP sexp_rdata2 = NULL;
	
	MPI_Comm_size(comm[commn], &gsize);
	displs=(int *)Calloc(gsize, int);
	displs[0]=0;
	for (i=1; i < gsize; i++)
		displs[i]=displs[i-1]+INTEGER(sexp_recvcounts)[i-1];

	switch (INTEGER(sexp_type)[0]){
	case 1:
		len=LENGTH(sexp_sdata);
		mpi_errhandler(MPI_Allgatherv(INTEGER(sexp_sdata), len, MPI_INT,
			INTEGER(sexp_rdata), INTEGER(sexp_recvcounts), displs,
				 MPI_INT,comm[commn]));
		break;
	case 2:
		len=LENGTH(sexp_sdata);
		mpi_errhandler(MPI_Allgatherv(REAL(sexp_sdata), len, MPI_DOUBLE,
			REAL(sexp_rdata), INTEGER(sexp_recvcounts), displs, 
				MPI_DOUBLE, comm[commn]));
		break;
	case 3:
	 	len=LENGTH(STRING_ELT(sexp_sdata,0));
                rlen=LENGTH(STRING_ELT(sexp_rdata,0));

                PROTECT (sexp_rdata2  = allocVector (STRSXP, 1));
                rdata = (char *)Calloc(rlen, char);
                MPI_Allgatherv(CHAR2 (STRING_ELT ((sexp_sdata),0)),len, MPI_CHAR, rdata,
		      INTEGER(sexp_recvcounts), displs, MPI_CHAR, comm[commn]);
                SET_STRING_ELT(sexp_rdata2, 0, mkChar(rdata));
                UNPROTECT(1);
		Free(rdata);
		break;
 	case 4:
                len=LENGTH(sexp_sdata);
                mpi_errhandler(MPI_Allgatherv(RAW(sexp_sdata), len, MPI_BYTE,
                        RAW(sexp_rdata), INTEGER(sexp_recvcounts), displs,
                                MPI_BYTE, comm[commn]));
                break;

	default:
		PROTECT(sexp_sdata=AS_NUMERIC(sexp_sdata));
		mpi_errhandler(MPI_Bcast(REAL(sexp_sdata), 1, datatype[0], 0, comm[commn]));
		UNPROTECT(1);
		break;	
	}
	Free(displs);
        if (INTEGER(sexp_type)[0]==3)
                return sexp_rdata2;
        else
                return sexp_rdata;
}

SEXP mpi_bcast(SEXP sexp_data,
			   SEXP sexp_type,
			   SEXP	sexp_rank,
			   SEXP sexp_comm,
			   SEXP sexp_buffunit){

	int len=LENGTH(sexp_data), type=INTEGER(sexp_type)[0];
	int rank=INTEGER(sexp_rank)[0], root,  commn=INTEGER(sexp_comm)[0],slen;
	int buffunit=INTEGER(sexp_buffunit)[0],errcode=0;
	char *rdata;
	SEXP sexp_data2 = NULL;
	//MPI_Datatype xdouble;
	R_xlen_t xlen=XLENGTH(sexp_data);
	
	switch (type){
	case 1:
		errcode=MPI_Bcast(INTEGER(sexp_data), len, MPI_INT, rank, comm[commn]);
		break;
	case 2:
		mpi_errhandler(MPI_Bcast(REAL(sexp_data), len, MPI_DOUBLE, rank, comm[commn]));
		break;
	case 3:
        	MPI_Comm_rank(comm[commn], &root);
		slen=LENGTH(STRING_ELT (sexp_data,0)); 
		if (rank==root) 
			MPI_Bcast(CHAR2 (STRING_ELT (sexp_data,0)), slen, 
				MPI_CHAR, rank, comm[commn]);
		else {
                	PROTECT (sexp_data2  = allocVector (STRSXP, 1));
	               	rdata = (char *)Calloc(slen, char);
                       	MPI_Bcast(rdata, slen, MPI_CHAR, rank, comm[commn]);
			SET_STRING_ELT(sexp_data2, 0, mkChar(rdata));
			UNPROTECT(1);
			Free(rdata);
		}
		break;
	case 4:
                errcode=MPI_Bcast(RAW(sexp_data), len, MPI_BYTE, rank, comm[commn]);
                break;
	case 5:
		MPI_Type_contiguous(buffunit, MPI_DOUBLE, xdouble);
		MPI_Type_commit(xdouble);
		if ((xlen % buffunit) > 0) len=1+(xlen/buffunit); else len=xlen/buffunit;
        mpi_errhandler(MPI_Bcast(REAL(sexp_data), len, xdouble[0], rank, comm[commn]));
		MPI_Type_free(xdouble);
		break;
	default:
		PROTECT(sexp_data=AS_NUMERIC(sexp_data));
		mpi_errhandler(MPI_Bcast(REAL(sexp_data), 1, datatype[0], rank, comm[commn]));
		UNPROTECT(1);
		break;		
	}
	if (errcode!=MPI_SUCCESS){
		int errmsglen;
		char errmsg[MPI_MAX_ERROR_STRING];
		MPI_Error_string(errcode, errmsg, &errmsglen);
		Rprintf("%s\n",errmsg);
		return mkString("error");
	}
	else {
        	if ((INTEGER(sexp_type)[0]==3) && (rank!=root))
                	return sexp_data2;
        	else
                	return sexp_data;
	}
}

SEXP mpi_send(SEXP sexp_data, 
			  SEXP sexp_type,
			  SEXP sexp_dest, 
			  SEXP sexp_tag,
			  SEXP sexp_comm){
	int slen,len=LENGTH(sexp_data),type=INTEGER(sexp_type)[0], dest=INTEGER(sexp_dest)[0];
	int commn=INTEGER(sexp_comm)[0], tag=INTEGER(sexp_tag)[0];

	switch (type){
	case 1:
		mpi_errhandler(MPI_Send(INTEGER(sexp_data), len, MPI_INT, dest, tag, comm[commn]));
		break;
	case 2:
		mpi_errhandler(MPI_Send(REAL(sexp_data), len, MPI_DOUBLE, dest, tag, comm[commn]));
		break;
	case 3:
		slen=LENGTH(STRING_ELT(sexp_data,0));
		MPI_Send(CHAR2(STRING_ELT(sexp_data,0)),slen, MPI_CHAR, dest, tag, comm[commn]); 
		break;
        case 4:
                MPI_Send(RAW(sexp_data),len, MPI_BYTE, dest, tag, comm[commn]);                
                break;

	default:
		PROTECT(sexp_data=AS_NUMERIC(sexp_data));
		mpi_errhandler(MPI_Send(REAL(sexp_data), 1, datatype[0], dest, tag, comm[commn]));
		UNPROTECT(1);
		break;		
	}
	return R_NilValue;
}

SEXP mpi_recv(SEXP sexp_data, 
  			  SEXP sexp_type,
			  SEXP sexp_source, 
			  SEXP sexp_tag,
			  SEXP sexp_comm,
			  SEXP sexp_status){
	int len=LENGTH(sexp_data), type=INTEGER(sexp_type)[0], source=INTEGER(sexp_source)[0];
	int tag=INTEGER(sexp_tag)[0],commn=INTEGER(sexp_comm)[0], statusn=INTEGER(sexp_status)[0];
	int slen;
	char *rdata;
	SEXP sexp_data2 = NULL;

	switch (type){
	case 1:
		mpi_errhandler(MPI_Recv(INTEGER(sexp_data), len, MPI_INT, source, tag, comm[commn],
			&status[statusn]));
		break;
	case 2:
		mpi_errhandler(MPI_Recv(REAL(sexp_data), len, MPI_DOUBLE, source, tag, comm[commn],
			&status[statusn]));
		break;
	case 3:
		slen=LENGTH(STRING_ELT(sexp_data,0));
                PROTECT (sexp_data2  = allocVector (STRSXP, 1));
                rdata = (char *)Calloc(slen, char);
		MPI_Recv(rdata, slen,MPI_CHAR,source,tag, comm[commn],&status[statusn]);
                SET_STRING_ELT(sexp_data2, 0, mkChar(rdata));
                UNPROTECT(1);
		Free(rdata);
		break;
	case 4:          
		mpi_errhandler(MPI_Recv(RAW(sexp_data), len, MPI_BYTE, source, tag, comm[commn],
                        &status[statusn]));
                break;

	default:
		PROTECT(sexp_data=AS_NUMERIC(sexp_data));
		mpi_errhandler(MPI_Recv(REAL(sexp_data), 1, datatype[0], source, tag, comm[commn],
			&status[statusn]));
		UNPROTECT(1);
		break;		
	}
        if (INTEGER(sexp_type)[0]==3)
                return sexp_data2;
        else
		return sexp_data;
}

SEXP mpi_reduce(SEXP sexp_send, 
				SEXP sexp_type,
				SEXP sexp_op, 
				SEXP sexp_dest,
				SEXP sexp_comm){
	int len=LENGTH(sexp_send), type=INTEGER(sexp_type)[0], dest=INTEGER(sexp_dest)[0];
	int commn=INTEGER(sexp_comm)[0], intop = INTEGER(sexp_op)[0];
	MPI_Op op= MPI_SUM;
	SEXP sexp_recv = NULL;

	switch(intop){
	case 1:
		op=MPI_SUM;
		break;
	case 2:
		op=MPI_PROD;
		break;
	case 3:
		op=MPI_MAX;
		break;
	case 4:
		op=MPI_MIN;
		break;
	case 5:
		op=MPI_MAXLOC;
		break;
	case 6:
		op=MPI_MINLOC;
		break;

	}
	switch(type){
	case 1:
		if (intop < 5){
			PROTECT (sexp_recv = allocVector(INTSXP, len));
			mpi_errhandler(MPI_Reduce(INTEGER(sexp_send), INTEGER(sexp_recv), 
			len, MPI_INT, op, dest, comm[commn])); 
			break;
		}
		else{
			int *send, rank, i;
			MPI_Comm_rank(comm[commn], &rank);
			send = (int *)Calloc(2*len, int);
			for (i=0; i < len; i++){
				send[2*i] = INTEGER(sexp_send)[i];
				send[2*i+1] = rank; 
			}
			PROTECT (sexp_recv = allocVector(INTSXP, 2*len));
			mpi_errhandler(MPI_Reduce(send, INTEGER(sexp_recv), 
			len, MPI_2INT, op, dest, comm[commn])); 
			Free(send);
			break;
		}
	case 2:
		if (intop < 5){
			PROTECT (sexp_recv = allocVector(REALSXP, len));
			mpi_errhandler(MPI_Reduce(REAL(sexp_send), REAL(sexp_recv), 
			len, MPI_DOUBLE, op, dest, comm[commn])); 
			break;
		}
		else {
			int i, rank;
			struct Dblint *send, *recv;
			send=(struct Dblint *)Calloc(len, struct Dblint);
			recv=(struct Dblint *)Calloc(len, struct Dblint);
			MPI_Comm_rank(comm[commn], &rank);
			for (i=0;i<len;i++){
				send[i].x = REAL(sexp_send)[i];
				send[i].rank = rank;
			}
			mpi_errhandler(MPI_Reduce(send, recv, len, MPI_DOUBLE_INT, op, dest, comm[commn])); 
			PROTECT (sexp_recv = allocVector(REALSXP, 2*len));
			for (i=0; i<len; i++){
				REAL(sexp_recv)[2*i] = recv[i].x;
				REAL(sexp_recv)[2*i+1] = recv[i].rank;
			}
			Free(send);
			Free(recv);
			break;
		}
	}
	
	UNPROTECT(1);			
	return sexp_recv;
}


SEXP mpi_allreduce(SEXP sexp_send, 
				   SEXP sexp_type,
				   SEXP sexp_op,
				   SEXP sexp_comm){
	int len=LENGTH(sexp_send), type=INTEGER(sexp_type)[0], commn=INTEGER(sexp_comm)[0];
	int intop = INTEGER(sexp_op)[0];
	MPI_Op op = MPI_SUM;
	SEXP sexp_recv = NULL;

	switch(intop){
	case 1:
		op=MPI_SUM;
		break;
	case 2:
		op=MPI_PROD;
		break;
	case 3:
		op=MPI_MAX;
		break;
	case 4:
		op=MPI_MIN;
		break;
	case 5:
		op=MPI_MAXLOC;
		break;
	case 6:
		op=MPI_MINLOC;
		break;

	}

	switch(type){
	case 1:
		if (intop < 5){
			PROTECT (sexp_recv = allocVector(INTSXP, len));
			mpi_errhandler(MPI_Allreduce(INTEGER(sexp_send), INTEGER(sexp_recv), 
			len, MPI_INT, op, comm[commn])); 
		break;
		}
		else{
			int *send, rank, i;
			MPI_Comm_rank(comm[commn], &rank);
			send = (int *)Calloc(2*len, int);
			for (i=0; i < len; i++){
				send[2*i] = INTEGER(sexp_send)[i];
				send[2*i+1] = rank; 
			}
			PROTECT (sexp_recv = allocVector(INTSXP, 2*len));
			mpi_errhandler(MPI_Allreduce(send, INTEGER(sexp_recv), 
			len, MPI_2INT, op, comm[commn])); 
			Free(send);
		break;
		}
	case 2:
		if (intop < 5) {
			PROTECT (sexp_recv = allocVector(REALSXP, len));
			mpi_errhandler(MPI_Allreduce(REAL(sexp_send), REAL(sexp_recv), 
				len, MPI_DOUBLE, op, comm[commn])); 
			break;
		}
		else {
			int i, rank;
			struct Dblint *send, *recv;
			send=(struct Dblint *)Calloc(len, struct Dblint);
			recv=(struct Dblint *)Calloc(len, struct Dblint);
			MPI_Comm_rank(comm[commn], &rank);
			for (i=0;i<len;i++){
				send[i].x = REAL(sexp_send)[i];
				send[i].rank = rank;
			}
			mpi_errhandler(MPI_Allreduce(send, recv, len, MPI_DOUBLE_INT, op, comm[commn])); 
			PROTECT (sexp_recv = allocVector(REALSXP, 2*len));
			for (i=0; i<len; i++){
				REAL(sexp_recv)[2*i] = recv[i].x;
				REAL(sexp_recv)[2*i+1] = recv[i].rank;
			}
			Free(send);
			Free(recv);
			break;
		}
	}
	
	UNPROTECT(1);			
	return sexp_recv;
}

SEXP mpi_iprobe(SEXP sexp_source, SEXP sexp_tag, SEXP sexp_comm, SEXP sexp_status){
	int flag;
	mpi_errhandler(MPI_Iprobe(INTEGER (sexp_source)[0], 
		INTEGER(sexp_tag)[0], comm[INTEGER(sexp_comm)[0]], &flag, 
		&status[INTEGER(sexp_status)[0]]));
	return AsInt(flag);
}

SEXP mpi_probe(SEXP sexp_source, SEXP sexp_tag, SEXP sexp_comm, SEXP sexp_status){
	return AsInt(erreturn(mpi_errhandler(MPI_Probe(INTEGER (sexp_source)[0], 
		INTEGER(sexp_tag)[0], comm[INTEGER(sexp_comm)[0]], 
		&status[INTEGER(sexp_status)[0]]))));
}

SEXP mpi_get_count(SEXP sexp_status, SEXP sexp_type){
	SEXP sexp_count;
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
	
	PROTECT (sexp_count = allocVector(INTSXP, 1));
	mpi_errhandler(MPI_Get_count(&status[INTEGER(sexp_status)[0]], datatype, INTEGER(sexp_count)));
	UNPROTECT(1);

	return sexp_count;
}

SEXP mpi_get_sourcetag (SEXP sexp_status){
	int statusn =INTEGER(sexp_status)[0];
	SEXP sexp_st;
	PROTECT(sexp_st=allocVector(INTSXP,2));
	INTEGER(sexp_st)[0]=status[statusn].MPI_SOURCE;
	INTEGER(sexp_st)[1]=status[statusn].MPI_TAG;
	UNPROTECT(1);
	return sexp_st;
}

/******************************* COMM **************************************/
SEXP mpi_barrier(SEXP sexp_comm){
	return AsInt(erreturn(mpi_errhandler(MPI_Barrier(comm[INTEGER(sexp_comm)[0]])))); 
}

SEXP mpi_comm_is_null(SEXP sexp_comm){
	return AsInt(comm[INTEGER(sexp_comm)[0]]==MPI_COMM_NULL);
}

SEXP mpi_comm_size(SEXP sexp_comm){
	int size;
	MPI_Comm_size(comm[INTEGER(sexp_comm)[0]], &size); 
	return AsInt(size);
}

SEXP mpi_comm_rank(SEXP sexp_comm){
	int rank;
	MPI_Comm_rank(comm[INTEGER(sexp_comm)[0]], &rank);
	return AsInt(rank);
}

SEXP mpi_comm_dup(SEXP sexp_comm, SEXP sexp_newcomm){
    int commn=INTEGER(sexp_comm)[0], newcommn=INTEGER(sexp_newcomm)[0];
    if (commn==0)
        return AsInt(erreturn(mpi_errhandler(MPI_Comm_dup(MPI_COMM_WORLD,
                &comm[newcommn]))));
    else
        return AsInt(erreturn(mpi_errhandler(MPI_Comm_dup(comm[commn],
                &comm[newcommn]))));
}

SEXP mpi_comm_c2f(SEXP sexp_comm){
  int c = INTEGER(sexp_comm)[0];
  return AsInt(MPI_Comm_c2f(comm[c]));
}

SEXP mpi_comm_free(SEXP sexp_comm){
	return AsInt(erreturn(mpi_errhandler(MPI_Comm_free(&comm[INTEGER(sexp_comm)[0]]))));
}

SEXP mpi_abort(SEXP sexp_comm){
	int errcode=0, commn=INTEGER(sexp_comm)[0];
	if (commn==0)
		MPI_Abort(MPI_COMM_WORLD, errcode);
	else
		MPI_Abort(comm[commn], errcode);
	Rprintf("The return errcode for mpi.abort() is %d\n", errcode);
	return AsInt(errcode);
}

/********************Intercomm********************************************/
SEXP mpi_comm_set_errhandler(SEXP sexp_comm){
	return AsInt(erreturn(MPI_Comm_set_errhandler(comm[INTEGER(sexp_comm)[0]], 
		MPI_ERRORS_RETURN)));
}

SEXP mpi_comm_test_inter(SEXP sexp_comm){
	int flag;
	MPI_Comm_test_inter(comm[INTEGER(sexp_comm)[0]], &flag);
	return AsInt(flag);
}

#ifdef MPI2
SEXP mpi_comm_spawn (SEXP sexp_slave,
					 SEXP sexp_argv,
					 SEXP sexp_nslave,
					 SEXP sexp_info,
					 SEXP sexp_root,
					 SEXP sexp_intercomm,
					 SEXP sexp_quiet){
    int i, nslave = INTEGER (sexp_nslave)[0], len = LENGTH (sexp_argv);
	int infon=INTEGER(sexp_info)[0], root=INTEGER(sexp_root)[0];
	int intercommn=INTEGER(sexp_intercomm)[0], *slaverrcode, realns;
    int quiet = INTEGER(sexp_quiet)[0];

	slaverrcode = (int *)Calloc(nslave, int);
	if (len==0)
		mpi_errhandler(MPI_Comm_spawn (CHAR2 (STRING_ELT (sexp_slave, 0)), MPI_ARGV_NULL, nslave,   
					info[infon], root, MPI_COMM_SELF, &comm[intercommn],
					slaverrcode)); 
	else {
		char **argv = (char **) R_alloc (len+1, sizeof (char *));
		for (i = 0; i < len; i++)
			argv[i] = CHAR2 (STRING_ELT (sexp_argv, i));
		argv[len] = NULL;
		mpi_errhandler(MPI_Comm_spawn (CHAR2 (STRING_ELT (sexp_slave, 0)), argv, nslave,   
					info[infon], root, MPI_COMM_SELF, &comm[intercommn],
					slaverrcode)); 
	}

	MPI_Comm_remote_size(comm[intercommn], &realns);
	if (realns < nslave)
		for (i=0; i < nslave; mpi_errhandler(slaverrcode[i++]));

	Free(slaverrcode);
	if (!quiet || realns < nslave)
		Rprintf("\t%d slaves are spawned successfully. %d failed.\n", realns, nslave-realns);
    return AsInt(realns);
}

SEXP mpi_comm_get_parent(SEXP sexp_comm){
	return AsInt(erreturn(mpi_errhandler(MPI_Comm_get_parent(&comm[INTEGER(sexp_comm)[0]]))));
}

SEXP mpi_is_master(){
	int check;
	MPI_Comm master;
	MPI_Comm_get_parent(&master);
	check=(master==MPI_COMM_NULL);
	if (!check) MPI_Comm_free(&master);
	return AsInt(check);
}

SEXP mpi_comm_disconnect(SEXP sexp_comm){
	return AsInt(erreturn(mpi_errhandler(MPI_Comm_disconnect(&comm[INTEGER(sexp_comm)[0]]))));
}
#endif

SEXP mpi_intercomm_merge(SEXP sexp_intercomm, SEXP sexp_high, SEXP sexp_comm){
	return AsInt(erreturn(mpi_errhandler(MPI_Intercomm_merge(comm[INTEGER(sexp_intercomm)[0]],
		INTEGER(sexp_high)[0],
		&comm[INTEGER(sexp_comm)[0]]))));
}

SEXP mpi_comm_remote_size(SEXP sexp_comm){
	int size;
	mpi_errhandler(MPI_Comm_remote_size(comm[INTEGER(sexp_comm)[0]], &size));
	return AsInt(size);
}

SEXP mpi_sendrecv(SEXP sexp_senddata,
        SEXP sexp_sendtype,
        SEXP sexp_dest,
        SEXP sexp_sendtag,
        SEXP sexp_recvdata,
        SEXP sexp_recvtype,                   
        SEXP sexp_source,                   
        SEXP sexp_recvtag,
        SEXP sexp_comm,
        SEXP sexp_status)
{
    int slen, rlen;
    int sendcount=LENGTH(sexp_senddata), sendtype=INTEGER(sexp_sendtype)[0];
    int dest=INTEGER(sexp_dest)[0], sendtag=INTEGER(sexp_sendtag)[0];
    int recvcount=LENGTH(sexp_recvdata), recvtype=INTEGER(sexp_recvtype)[0];
    int source=INTEGER(sexp_source)[0], recvtag=INTEGER(sexp_recvtag)[0];
    int commn=INTEGER(sexp_comm)[0],statusn=INTEGER(sexp_status)[0];
    char *rdata;
    SEXP sexp_recvdata2 = NULL;

    switch(sendtype){
        case 1:
            switch(recvtype){
               	case 1:
                    MPI_Sendrecv(INTEGER(sexp_senddata),sendcount, 
			MPI_INT,dest,sendtag, INTEGER(sexp_recvdata), 
			recvcount, MPI_INT, source,recvtag,
                        comm[commn], &status[statusn]);
                    break;
                case 2:
                    MPI_Sendrecv(INTEGER(sexp_senddata),sendcount, 
			MPI_INT,dest,sendtag, REAL(sexp_recvdata), 
			recvcount, MPI_DOUBLE, source,recvtag,
                        comm[commn], &status[statusn]);
                    break;
                case 3:
                    rlen=LENGTH(STRING_ELT(sexp_recvdata,0)); 
                    PROTECT (sexp_recvdata2  = allocVector (STRSXP, 1));
                    rdata = (char *)Calloc(rlen, char);
                    MPI_Sendrecv(INTEGER(sexp_senddata),sendcount, MPI_INT, dest, sendtag, 
			rdata, rlen, MPI_CHAR, source, recvtag, comm[commn],  &status[statusn]);
                    SET_STRING_ELT(sexp_recvdata2, 0, mkChar(rdata));
                    UNPROTECT(1);
                    Free(rdata);
                    break;
 		case 4:
                     MPI_Sendrecv(INTEGER(sexp_senddata),sendcount,
                        MPI_INT,dest,sendtag, RAW(sexp_recvdata),
                        recvcount, MPI_BYTE, source,recvtag,
                        comm[commn], &status[statusn]);
		    break;
                }
	    break;
        case 2:
            switch(recvtype){
                case 1:
                    MPI_Sendrecv(REAL(sexp_senddata),sendcount, 
			MPI_DOUBLE,dest,sendtag, INTEGER(sexp_recvdata), 
			recvcount, MPI_INT, source,recvtag,
                        comm[commn], &status[statusn]);
                    break;
                case 2: 
                    MPI_Sendrecv(REAL(sexp_senddata),sendcount, 
			MPI_DOUBLE,dest,sendtag, REAL(sexp_recvdata), 
			recvcount, MPI_DOUBLE, source,recvtag,
                        comm[commn], &status[statusn]);
                    break;
                case 3:
                    rlen=LENGTH(STRING_ELT(sexp_recvdata,0)); 
                    PROTECT (sexp_recvdata2  = allocVector (STRSXP, 1));
                    rdata = (char *)Calloc(rlen, char);
                    MPI_Sendrecv(REAL(sexp_senddata),sendcount, MPI_DOUBLE, dest, sendtag, 
			rdata, rlen, MPI_CHAR, source, recvtag, comm[commn],  &status[statusn]);
                    SET_STRING_ELT(sexp_recvdata2, 0, mkChar(rdata));
                    UNPROTECT(1);
                    Free(rdata);
                    break;
 		case 4:
                    MPI_Sendrecv(REAL(sexp_senddata),sendcount,
                        MPI_DOUBLE,dest,sendtag, RAW(sexp_recvdata),
                        recvcount, MPI_BYTE, source,recvtag,
                        comm[commn], &status[statusn]);
                    break;

                }
            break;
        case 3: 
            slen=LENGTH(STRING_ELT(sexp_senddata,0));                       
  
            switch(recvtype){
               	case 1:
                    MPI_Sendrecv(CHAR2(STRING_ELT(sexp_senddata,0)),slen, 
			MPI_CHAR,dest,sendtag, INTEGER(sexp_recvdata), 
			recvcount, MPI_INT, source,recvtag,
                        comm[commn], &status[statusn]);
                    break;
                case 2:                       
		    MPI_Sendrecv(CHAR2(STRING_ELT(sexp_senddata,0)),slen, 
			MPI_CHAR,dest,sendtag, REAL(sexp_recvdata), 
			recvcount, MPI_DOUBLE, source,recvtag,
                        comm[commn], &status[statusn]);
                    break;
                case 3:
                    rlen=LENGTH(STRING_ELT(sexp_recvdata,0));
                    PROTECT (sexp_recvdata2  = allocVector (STRSXP, 1));
                    rdata = (char *)Calloc(rlen, char);

		    MPI_Sendrecv(CHAR2(STRING_ELT(sexp_senddata,0)),slen, MPI_CHAR, dest, sendtag, 
			rdata, rlen, MPI_CHAR, source, recvtag, comm[commn], &status[statusn]);
                    SET_STRING_ELT(sexp_recvdata2, 0, mkChar(rdata));
                    UNPROTECT(1);
                    Free(rdata);
                    break;
 		case 4:
                    MPI_Sendrecv(CHAR2(STRING_ELT(sexp_senddata,0)),slen,
                        MPI_CHAR,dest,sendtag, RAW(sexp_recvdata),
                        recvcount, MPI_BYTE, source,recvtag,
                        comm[commn], &status[statusn]);
                    break;

          	}
   		break;
 	case 4:
            switch(recvtype){
                case 1:
                    MPI_Sendrecv(RAW(sexp_senddata),sendcount,
                        MPI_BYTE,dest,sendtag, INTEGER(sexp_recvdata),
                        recvcount, MPI_INT, source,recvtag,
                        comm[commn], &status[statusn]);
                    break;
                case 2:
                    MPI_Sendrecv(RAW(sexp_senddata),sendcount,
                        MPI_BYTE,dest,sendtag, REAL(sexp_recvdata),
                        recvcount, MPI_DOUBLE, source,recvtag,
                        comm[commn], &status[statusn]);
                    break;
                case 3:
                    rlen=LENGTH(STRING_ELT(sexp_recvdata,0));
                    PROTECT (sexp_recvdata2  = allocVector (STRSXP, 1));
                    rdata = (char *)Calloc(rlen, char);
                    MPI_Sendrecv(RAW(sexp_senddata),sendcount, MPI_BYTE, dest, sendtag, 
			rdata, rlen, MPI_CHAR, source, recvtag, comm[commn],  &status[statusn]);
                    SET_STRING_ELT(sexp_recvdata2, 0, mkChar(rdata));
                    UNPROTECT(1);
                    Free(rdata);
                    break;
                case 4:
                    MPI_Sendrecv(RAW(sexp_senddata),sendcount,
                        MPI_BYTE,dest,sendtag, RAW(sexp_recvdata),
                        recvcount, MPI_BYTE, source,recvtag,
                        comm[commn], &status[statusn]);
                    break;
                }
            break;

	    }
    if (recvtype==3)
	return sexp_recvdata2;
    else
    	return sexp_recvdata;          
}

SEXP mpi_sendrecv_replace(SEXP sexp_data,
        SEXP sexp_type,
        SEXP sexp_dest,
        SEXP sexp_sendtag,
        SEXP sexp_source,
        SEXP sexp_recvtag,
        SEXP sexp_comm, 
        SEXP sexp_status)
{
        int slen;
        int len=LENGTH(sexp_data), type=INTEGER(sexp_type)[0];
        int dest=INTEGER(sexp_dest)[0], sendtag=INTEGER(sexp_sendtag)[0];
        int source=INTEGER(sexp_source)[0],recvtag=INTEGER(sexp_recvtag)[0];
        int commn=INTEGER(sexp_comm)[0],statusn=INTEGER(sexp_status)[0];
	char *srdata;
	SEXP sexp_data2 = NULL;

        switch (type){
        case 1:
                MPI_Sendrecv_replace(INTEGER(sexp_data), len, MPI_INT, dest, 
		sendtag, source, recvtag, comm[commn], &status[statusn]);
                break;
        case 2:
                MPI_Sendrecv_replace(REAL(sexp_data), len, MPI_DOUBLE, dest, 
		sendtag, source, recvtag, comm[commn], &status[statusn]);
                break;

        case 3:
                slen=LENGTH(STRING_ELT(sexp_data,0));
		PROTECT (sexp_data2  = allocVector (STRSXP, 1));
		srdata= (char *)Calloc(slen, char);
		strcpy(srdata, CHAR(STRING_ELT(sexp_data,0)));
                MPI_Sendrecv_replace(srdata, slen,MPI_CHAR, dest, sendtag, source, recvtag, 
				comm[commn], &status[statusn]); 
		UNPROTECT(1);
		Free(srdata); 
                break;
 	case 4:
                MPI_Sendrecv_replace(RAW(sexp_data), len, MPI_BYTE, dest,
                sendtag, source, recvtag, comm[commn], &status[statusn]);
                break;

        default:
                PROTECT(sexp_data=AS_NUMERIC(sexp_data));
                MPI_Sendrecv_replace(REAL(sexp_data), 1, datatype[0], dest, 
			sendtag, source, recvtag, comm[commn],
                        &status[statusn]);                        
                break;
          }
	if (type==3)
		return sexp_data2;
	else
          	return sexp_data;
}

/************ cart dim *************************************/

SEXP mpi_cart_create(SEXP sexp_comm_old,  SEXP sexp_dims, SEXP sexp_periods, SEXP sexp_reorder, 
           SEXP sexp_comm_cart) {
        int comm_old = INTEGER(sexp_comm_old)[0];
        int ndims = LENGTH(sexp_dims);
        int reorder = INTEGER(sexp_reorder)[0];
        int comm_cart = INTEGER(sexp_comm_cart)[0];
        int retcode; 
        retcode=erreturn(mpi_errhandler(MPI_Cart_create(comm[comm_old], ndims, 
                INTEGER(sexp_dims), INTEGER(sexp_periods), reorder, &comm[comm_cart])));    
        return  AsInt(retcode);
}

SEXP mpi_dims_create(SEXP sexp_nnodes, SEXP sexp_ndims, SEXP sexp_dims) {
        int nnodes = INTEGER(sexp_nnodes)[0];
        int ndims = INTEGER(sexp_ndims)[0];
        mpi_errhandler(MPI_Dims_create(nnodes, ndims, INTEGER(sexp_dims)));
        return sexp_dims;
}


SEXP mpi_cartdim_get(SEXP sexp_comm) {
        int comm2 = INTEGER(sexp_comm)[0];
        int ndims;
        mpi_errhandler(MPI_Cartdim_get(comm[comm2], &ndims));
        return AsInt(ndims);    
}

SEXP mpi_cart_get(SEXP sexp_comm, SEXP sexp_maxdims) {
        int comm2 = INTEGER(sexp_comm)[0];
        int maxdims = INTEGER(sexp_maxdims)[0];
        SEXP dims_periods_coords;

        PROTECT (dims_periods_coords = allocVector(INTSXP, maxdims*3));
        
        mpi_errhandler(MPI_Cart_get(comm[comm2], maxdims, INTEGER(dims_periods_coords), 
		INTEGER(dims_periods_coords) + maxdims, INTEGER(dims_periods_coords) + maxdims*2));
        
        UNPROTECT(1);
        return dims_periods_coords;
}


SEXP mpi_cart_rank(SEXP sexp_comm, SEXP sexp_coords){
        int comm2 = INTEGER(sexp_comm)[0];      
        int rank;
        mpi_errhandler(MPI_Cart_rank(comm[comm2], INTEGER(sexp_coords), &rank));
        return AsInt(rank);
}

SEXP mpi_cart_coords(SEXP sexp_comm, SEXP sexp_rank, SEXP sexp_maxdims) {
        int comm2 = INTEGER(sexp_comm)[0];
        int rank = INTEGER(sexp_rank)[0];
        int maxdims = INTEGER(sexp_maxdims)[0];
        SEXP coords;
        PROTECT (coords = allocVector(INTSXP, maxdims));
        mpi_errhandler(MPI_Cart_coords(comm[comm2], rank, maxdims, INTEGER(coords)));
        UNPROTECT(1);   
        return coords;
}


SEXP mpi_cart_shift(SEXP sexp_comm, SEXP sexp_direction, SEXP sexp_disp) {
        int comm2 = INTEGER(sexp_comm)[0];
        int direction = INTEGER(sexp_direction)[0];
        int disp = INTEGER(sexp_disp)[0];
        SEXP rank_source_dest;  
        PROTECT (rank_source_dest = allocVector(INTSXP,2 ));
        mpi_errhandler(MPI_Cart_shift(comm[comm2], direction, disp, &INTEGER(rank_source_dest)[0],
		 &INTEGER(rank_source_dest)[1]));
        UNPROTECT(1);
        return rank_source_dest;
}

/************** nonblocking point to point calls *************************/
SEXP mpi_isend(SEXP sexp_data,
                          SEXP sexp_type,
                          SEXP sexp_dest,
                          SEXP sexp_tag,
                          SEXP sexp_comm,
                          SEXP sexp_request){
        int slen,len=LENGTH(sexp_data),type=INTEGER(sexp_type)[0], dest=INTEGER(sexp_dest)[0];
        int commn=INTEGER(sexp_comm)[0], tag=INTEGER(sexp_tag)[0], requestn=INTEGER(sexp_request)[0];
        
        switch (type){
        case 1:
                mpi_errhandler(MPI_Isend(INTEGER(sexp_data), len, MPI_INT, dest, tag, 
			comm[commn], &request[requestn]));
                break;
        case 2:

                mpi_errhandler(MPI_Isend(REAL(sexp_data), len, MPI_DOUBLE, dest, tag, comm[commn], 
			&request[requestn]));
                break;
        case 3:
                slen=LENGTH(STRING_ELT(sexp_data,0));
                mpi_errhandler(MPI_Isend(CHAR2(STRING_ELT(sexp_data,0)),slen, MPI_CHAR, dest, 
			tag, comm[commn], &request[requestn]));
                break;
	case 4:

                mpi_errhandler(MPI_Isend(RAW(sexp_data), len, MPI_BYTE, dest, tag, comm[commn],
                        &request[requestn]));
                break;

        default:
                PROTECT(sexp_data=AS_NUMERIC(sexp_data));
                mpi_errhandler(MPI_Isend(REAL(sexp_data), 1, datatype[0], dest, tag, comm[commn], 
			&request[requestn]));
                UNPROTECT(1);
                break;
        }
        return R_NilValue;
}

SEXP mpi_irecv(SEXP sexp_data,
                          SEXP sexp_type,
                          SEXP sexp_source,
                          SEXP sexp_tag,
                          SEXP sexp_comm,
                          SEXP sexp_request){
        int slen,len=LENGTH(sexp_data),type=INTEGER(sexp_type)[0], source=INTEGER(sexp_source)[0];
        int commn=INTEGER(sexp_comm)[0], tag=INTEGER(sexp_tag)[0], requestn=INTEGER(sexp_request)[0];

        switch (type){
        case 1:
                mpi_errhandler(MPI_Irecv(INTEGER(sexp_data), len, MPI_INT, source, tag, 
			comm[commn], &request[requestn]));
                break;
        case 2:
                mpi_errhandler(MPI_Irecv(REAL(sexp_data), len, MPI_DOUBLE, source, tag, comm[commn], 
			&request[requestn]));
                break;
        case 3:
                slen=LENGTH(STRING_ELT(sexp_data,0));
                        mpi_errhandler(MPI_Irecv(CHAR2(STRING_ELT(sexp_data,0)),slen, MPI_CHAR, 
				source, tag, comm[commn], &request[requestn]));
                break;
 	case 4:
                mpi_errhandler(MPI_Irecv(RAW(sexp_data), len, MPI_BYTE, source, tag, comm[commn],
                        &request[requestn]));
                break;

        default:
                PROTECT(sexp_data=AS_NUMERIC(sexp_data));
                mpi_errhandler(MPI_Irecv(REAL(sexp_data), 1, datatype[0], source, tag, comm[commn], 
			&request[requestn]));
                UNPROTECT(1);
                break;
        }
        return R_NilValue;
}

SEXP mpi_wait(SEXP sexp_request, SEXP sexp_status){
        int requestn=INTEGER(sexp_request)[0], statusn=INTEGER(sexp_status)[0];
        mpi_errhandler(MPI_Wait(&request[requestn], &status[statusn]));
        return R_NilValue;
}


SEXP mpi_test(SEXP sexp_request,  SEXP sexp_status){
        int requestn=INTEGER(sexp_request)[0], flag, statusn=INTEGER(sexp_status)[0];
        mpi_errhandler(MPI_Test(&request[requestn], &flag, &status[statusn]));
        return AsInt(flag);
}

SEXP mpi_cancel(SEXP sexp_request){
	int requestn=INTEGER(sexp_request)[0];
	mpi_errhandler(MPI_Cancel(&request[requestn]));
	return R_NilValue;
}

SEXP mpi_test_cancelled(SEXP sexp_status){
        int flag, statusn=INTEGER(sexp_status)[0];
        mpi_errhandler(MPI_Test_cancelled(&status[statusn], &flag));
        return AsInt(flag);
}

SEXP mpi_waitany(SEXP sexp_count, SEXP sexp_status){
        int index, countn=INTEGER(sexp_count)[0],statusn=INTEGER(sexp_status)[0];
        mpi_errhandler(MPI_Waitany(countn, request, &index, &status[statusn]));
        return AsInt(index);
}

SEXP mpi_testany(SEXP sexp_count, SEXP sexp_status){
        int countn=INTEGER(sexp_count)[0],  statusn=INTEGER(sexp_status)[0];
		SEXP indexflag;
		PROTECT (indexflag = allocVector(INTSXP, 2));
        mpi_errhandler(MPI_Testany(countn, request, &INTEGER(indexflag)[0],
		&INTEGER(indexflag)[1], &status[statusn]));
		UNPROTECT(1);
        return indexflag;
}

SEXP mpi_waitall(SEXP sexp_count){
        int countn=INTEGER(sexp_count)[0];
        mpi_errhandler(MPI_Waitall(countn, request, status));
        return R_NilValue;
}

SEXP mpi_testall(SEXP sexp_count){
        int countn=INTEGER(sexp_count)[0], flag;
        mpi_errhandler(MPI_Testall(countn, request, &flag, status));
        return AsInt(flag);
}

SEXP mpi_testsome(SEXP sexp_count){
        int countn=INTEGER(sexp_count)[0];
		SEXP indices;
		PROTECT (indices = allocVector(INTSXP, countn+1));
        mpi_errhandler(MPI_Testsome(countn, request, &INTEGER(indices)[0], 
		&INTEGER(indices)[1], status));
		UNPROTECT(1);
        return indices;
}

SEXP mpi_waitsome(SEXP sexp_count){
        int countn=INTEGER(sexp_count)[0];
		SEXP indices;
		PROTECT (indices = allocVector(INTSXP, countn+1));
        mpi_errhandler(MPI_Waitsome(countn, request, &INTEGER(indices)[0], 
		&INTEGER(indices)[1], status));
		UNPROTECT(1);
        return indices;
}

/*
SEXP mpi_request_get_status(SEXP sexp_request,  SEXP sexp_status){
        int requestn=INTEGER(sexp_request)[0], flag, statusn=INTEGER(sexp_status)[0];
        mpi_errhandler(MPI_Request_get_status(request[requestn], &flag, &status[statusn]));
        return AsInt(flag);
}
*/

