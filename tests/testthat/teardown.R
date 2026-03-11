if (!identical(Sys.getenv("_R_CHECK_PACKAGE_NAME_", ""), "npRmpi")) {
  try(close_mpi_slaves(), silent = TRUE)
}
