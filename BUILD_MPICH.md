# npRmpi Build/Run (macOS + MacPorts + MPICH)

This is the working recipe for local build/check and runtime sanity checks.

## 1) Install MPICH with MacPorts

```bash
sudo port install mpich-default
sudo port select --set mpi mpich-mp-fortran
```

## 2) Build environment

```bash
export RMPI_TYPE=MPICH
export RMPI_INCLUDE=/opt/local/include/mpich-mp
export RMPI_LIB_PATH=/opt/local/lib/mpich-mp
export RMPI_LIBS="-L/opt/local/lib/mpich-mp -lmpi"
export CC=mpicc
export CXX=mpicxx
```

## 3) Build/install/check

```bash
cd /Users/jracine/Development
R CMD build np-npRmpi
R CMD INSTALL npRmpi_0.70-0.tar.gz
R CMD check --as-cran npRmpi_0.70-0.tar.gz
```

## 4) Interactive runtime smoke

```r
library(npRmpi)
npRmpi.start(mode = "spawn", nslaves = 1)
options(npRmpi.autodispatch = TRUE, np.messages = FALSE)

set.seed(1)
x <- runif(200)
y <- sin(2*pi*x) + rnorm(200, sd = 0.2)
bw <- npregbw(y ~ x, regtype = "ll", bwmethod = "cv.ls")
fit <- npreg(bws = bw)
summary(fit)

npRmpi.stop()
```

## 5) Batch/cluster runtime (no external profile bootstrap)

`foo.R` should include:

```r
library(npRmpi)
npRmpi.start(mode = "attach", autodispatch = TRUE, np.messages = FALSE)
# ... np* calls ...
npRmpi.stop(mode = "attach")
```

Then launch:

```bash
mpiexec -n 128 Rscript foo.R
```

## Notes

- If `MPI_Init` fails during build/check in restricted environments, set:

```bash
export NP_RMPI_SKIP_INIT=1
```

- For sandboxed MPI runs on macOS, use:

```bash
FI_TCP_IFACE=en0
```
