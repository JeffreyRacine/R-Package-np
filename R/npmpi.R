np.mpi.initialize = function(){
  if ((mpi.status=.C("np_mpi_init", mpi.status = integer(2), PACKAGE= "npRmpi")$mpi.status)[1] == -1)
    stop("MPI support has not been enabled. Please recompile np with --enable-MPI")
  invisible(mpi.status)
}
