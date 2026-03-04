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

### MPI Example Modes During Check

- Default `R CMD check --as-cran` behavior:
  - MPI spawn examples are skipped in check context to avoid MPICH teardown
    killing the parent check process on some systems.
- To force MPI spawn examples to run during check:

```bash
cd /Users/jracine/Development/np-npRmpi/man
./run

cd /Users/jracine/Development
NP_RMPI_RUN_MPI_EXAMPLES_IN_CHECK=1 R CMD check --as-cran npRmpi_0.70-0.tar.gz

cd /Users/jracine/Development/np-npRmpi/man
./dontrun
```

## Runtime Modes

- Interactive R session:
  - `npRmpi.init(mode="spawn", nslaves=...)`
- Cluster/batch under `mpiexec`:
  - start script with `npRmpi.init(mode="attach", autodispatch=TRUE, np.messages=FALSE)`
  - end script with `npRmpi.quit(mode="attach")`
