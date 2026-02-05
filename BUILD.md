# Build (npRmpi)

## Environment (MPICH via MacPorts)

```bash
export RMPI_TYPE=MPICH
export RMPI_INCLUDE=/opt/local/include/mpich-mp
export RMPI_LIB_PATH=/opt/local/lib/mpich-mp
export RMPI_LIBS="-L/opt/local/lib/mpich-mp -lmpi"
export CC=mpicc
export CXX=mpicxx
```

## Build Tarball

```bash
cd /Users/jracine/Development/np-npRmpi
R CMD build np
```

This produces `npRmpi_0.60-21.tar.gz` in the same directory.

## Install

```bash
cd /Users/jracine/Development/np-npRmpi
R CMD INSTALL npRmpi_0.60-21.tar.gz
```

## Quick Load Check

```bash
R -q -e 'library(npRmpi); sessionInfo()'
```
