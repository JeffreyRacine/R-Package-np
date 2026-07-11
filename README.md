# npRmpi

`npRmpi` provides parallel nonparametric kernel methods for mixed datatypes.
It is written and maintained by Jeffrey S. Racine
(racinej@mcmaster.ca) and co-authored by Tristen Hayfield
(tristen.hayfield@gmail.com).

## Installation

With a working MPI implementation already installed, install the current CRAN
release with:

```r
install.packages("npRmpi", dependencies = TRUE)
```

The CRAN release and the development branch can have different version
numbers. To install the development branch directly from GitHub:

```r
remotes::install_github(
  "JeffreyRacine/R-Package-np",
  ref = "npRmpi",
  build_vignettes = FALSE
)
```

Set `build_vignettes = TRUE` only when the required vignette toolchain is
available.

## Supported MPI backend for 0.70-5

The 0.70-5 release protocol is validated with MPICH. Open MPI runtime
validation is explicitly waived for this release because its ordinary
singleton-spawn and prelaunched-master finalization behavior did not satisfy
the release gates. Do not assume that an Open-MPI-linked build has the same
supported runtime contract.

### MPICH via MacPorts

```bash
sudo port install mpich-default
sudo port select --set mpi mpich-mp-fortran

export RMPI_TYPE=MPICH
export RMPI_INCLUDE=/opt/local/include/mpich-mp
export RMPI_LIB_PATH=/opt/local/lib/mpich-mp
```

The package configure script consumes these three `RMPI_*` variables.
`RMPI_LIBS` is not required. Let R use its configured C/C++/Fortran compilers;
do not replace `CC` and `CXX` with MPI wrapper compilers unless you have
deliberately validated that toolchain combination.

Verify the selected runtime before building:

```bash
mpiexec --version
test -f "$RMPI_INCLUDE/mpi.h"
test -f "$RMPI_LIB_PATH/libmpi.dylib"
R CMD config CC
R CMD config CXX
R CMD config FLIBS
```

From the repository root, build and install a source tarball without
hardcoding the package version:

```bash
VERSION=$(awk -F': *' '/^Version:/{print $2; exit}' DESCRIPTION)
R CMD build .
R CMD INSTALL "npRmpi_${VERSION}.tar.gz"
```

See [BUILD_MPICH.md](BUILD_MPICH.md) for the complete public MPICH recipe and
[BUILD.md](BUILD.md) for release-validation expectations.

## Quick start

Interactive session/spawn mode:

```r
library(npRmpi)
npRmpi.init(nslaves = 1)

set.seed(1)
x <- runif(200)
y <- sin(2 * pi * x) + rnorm(200, sd = 0.2)
bw <- npregbw(y ~ x, regtype = "ll", bwmethod = "cv.ls")
fit <- npreg(bws = bw)
summary(fit)

npRmpi.quit(force = TRUE)
```

On macOS, `npRmpi.quit()` without `force = TRUE` may intentionally retain the
worker pool for reuse. Use `force = TRUE` when the session is finished and a
clean teardown is required.

Batch/cluster attach mode:

```bash
mpiexec -n 4 Rscript foo.R
```

`foo.R` should contain:

```r
library(npRmpi)
is_master <- npRmpi.init(mode = "attach", autodispatch = TRUE,
                         np.messages = FALSE)

if (isTRUE(is_master)) {
  # ... np* calls ...
  npRmpi.quit(mode = "attach")
}
```

See the installed `MPI_SETUP.md` guide for session, attach, and
profile/manual-broadcast details.

For more information, visit the maintainer's website:
<https://experts.mcmaster.ca/people/racinej>.
