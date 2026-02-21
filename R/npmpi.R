np.mpi.initialize = function(){
  mpi.status <- .Call("C_np_mpi_init", PACKAGE= "npRmpi")
  if (mpi.status[1] == -1)
    stop("MPI support has not been enabled. Please recompile np with --enable-MPI")
  invisible(mpi.status)
}
