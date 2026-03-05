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
npRmpi.init(mode = "spawn", nslaves = 1)
options(npRmpi.autodispatch = TRUE, np.messages = FALSE)

set.seed(1)
x <- runif(200)
y <- sin(2*pi*x) + rnorm(200, sd = 0.2)
bw <- npregbw(y ~ x, regtype = "ll", bwmethod = "cv.ls")
fit <- npreg(bws = bw)
summary(fit)

npRmpi.quit()
```

## 5) Batch/cluster runtime (no external profile bootstrap)

`foo.R` should include:

```r
library(npRmpi)
npRmpi.init(mode = "attach", autodispatch = TRUE, np.messages = FALSE)
# ... np* calls ...
npRmpi.quit(mode = "attach")
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

- If `R CMD check` still fails in sandbox with OFI/NIC initialization errors, use the loopback + tcp provider fallback:

```bash
FI_TCP_IFACE=lo0 FI_PROVIDER=tcp FI_SOCKETS_IFACE=lo0 R CMD check --no-manual npRmpi_0.60-21.tar.gz
```

## Canonical Implementation Directive (2026-03-05)

This repository follows a strict canonical execution rule:

1. One canonical implementation per method (outside explicit `np.tree` branching).
2. Unsupported configurations must fail fast with explicit `stop(...)` diagnostics.
3. No silent remap/coercion of user-selected options (for example `bwmethod`, `regtype`, kernels, `cv.iterate`, or bounds transforms).
4. No hidden alternate execution paths for the same method semantics.
5. All fit-defining options (for example `degree`, `basis`, `bernstein.basis`, kernels, and bounds) must be propagated and used by the canonical path.
6. `np.tree=FALSE` is the default; when `np.tree=TRUE`, behavior must remain semantics-preserving and option-compatible with the canonical path.
7. Remove or reject legacy/debug compatibility branches that add redundant runtime overhead once canonical behavior exists.
