#include <mpi.h>
#include <stdio.h>
#include <string.h>

int main(int argc, char **argv){
        int i, *univ_size, univ_flag, *n;
        int rank, size;
	double mypi, pi, h, sum, x;

        MPI_Comm slavecomm, all_processes;

        /*Initialize MPI*/
        MPI_Init(&argc, &argv); 
        
        MPI_Comm_get_parent(&slavecomm);
        MPI_Intercomm_merge(slavecomm, 1, &all_processes);

        /*How many processes are there?*/
        MPI_Comm_size(all_processes, &size);
 
        /*Which one am I?*/
        MPI_Comm_rank(all_processes, &rank);

        MPI_Bcast(n,1,MPI_INT,0, all_processes);

        /*Compute portion of pi on each node */
        h   = 1.0 / (double) (*n);
        sum = 0.0;
        for (i = rank ; i <= *n; i += size-1) {
         x = h * ((double)i - 0.5);
            sum += (4.0 / (1.0 + x*x));
         }
        mypi = h * sum;

	/*Send all computed portion of pi to rank 0 and sum together */
	MPI_Reduce(&mypi, &pi,1, MPI_DOUBLE, MPI_SUM,0, all_processes);

        /*All done*/
	MPI_Comm_free(&all_processes);
        MPI_Finalize();
        exit(0);
}	
