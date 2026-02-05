# npRmpi build/run (macOS + MacPorts + MPICH)

This is the working recipe to build and run `npRmpi` on macOS using MacPorts + MPICH.

## 1) Install MPICH with MacPorts
```bash
sudo port install mpich-default
sudo port select --set mpi mpich-mp-fortran
```

## 2) Environment variables (required for build)
```bash
export RMPI_TYPE=MPICH
export RMPI_INCLUDE=/opt/local/include/mpich-mp
export RMPI_LIB_PATH=/opt/local/lib/mpich-mp
export RMPI_LIBS="-L/opt/local/lib/mpich-mp -lmpi"
export CC=mpicc
export CXX=mpicxx
```

## 3) Build and install npRmpi
```bash
cd ~/Development
R CMD build np
R CMD INSTALL npRmpi_0.60-21.tar.gz
```

## 4) Quick runtime sanity check
```r
library(npRmpi)
library(Rmpi)
mpi.spawn.Rslaves(nslaves=1)
mpi.remote.exec(paste("Hello from", mpi.comm.rank()))
mpi.close.Rslaves()
mpi.quit()
```

## Notes
- If you see an MPI_Init error during `R CMD build` (e.g. OFI/`en1` permissions), you can set:
  ```bash
  export NP_RMPI_SKIP_INIT=1
  ```
  This skips MPI initialization during package load for the build stage only.
- If you previously used OpenMPI, ensure MPICH is active via `port select --set mpi mpich-mp-fortran`.
- The build expects MPI2 support. The configure script now adds `-DMPI2` for MPICH and links against `-lmpi`.
