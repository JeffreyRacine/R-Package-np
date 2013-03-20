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


int mpi_errhandler(int errcode)
{
	int errmsglen;
	char errmsg[MPI_MAX_ERROR_STRING];

    if (errcode != MPI_SUCCESS) {
		MPI_Error_string(errcode, errmsg, &errmsglen);
		error(errmsg);
	}
	
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
