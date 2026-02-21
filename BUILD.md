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
cd /Users/jracine/Development
R CMD build np-npRmpi
```

This produces `npRmpi_0.70-0.tar.gz` in `/Users/jracine/Development`.

## Install

```bash
cd /Users/jracine/Development
R CMD INSTALL npRmpi_0.70-0.tar.gz
```

## Quick Load Check

```bash
R -q -e 'library(npRmpi); sessionInfo()'
```

## Check

```bash
cd /Users/jracine/Development
R CMD check --as-cran npRmpi_0.70-0.tar.gz
```

## Runtime Modes

- Interactive R session:
  - `npRmpi.init(mode="spawn", nslaves=...)`
- Cluster/batch under `mpiexec`:
  - start script with `npRmpi.init(mode="attach", autodispatch=TRUE, np.messages=FALSE)`
  - end script with `npRmpi.quit(mode="attach")`
