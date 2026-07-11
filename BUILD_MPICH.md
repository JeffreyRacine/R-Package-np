# Building npRmpi with MPICH on macOS

This is the authoritative public build recipe for the supported npRmpi 0.70-5
backend on macOS. It was validated with MacPorts MPICH 4.3.2 on Apple Silicon.
A future MPICH major-version upgrade requires fresh installed-tarball runtime,
parity, lifecycle, and scaling validation.

Open MPI is not a supported release backend for 0.70-5; its runtime validation
is explicitly waived.

## 1. Install and select MPICH

```bash
sudo port install mpich-default
sudo port select --set mpi mpich-mp-fortran
```

Confirm that MacPorts selected MPICH rather than another MPI implementation:

```bash
port select --list mpi
command -v mpicc
command -v mpiexec
mpiexec --version
```

The selected `mpicc`, `mpicxx`, and `mpiexec` commands should resolve through
MacPorts to the `*-mpich-mp` wrappers.

## 2. Configure MPI discovery

```bash
export RMPI_TYPE=MPICH
export RMPI_INCLUDE=/opt/local/include/mpich-mp
export RMPI_LIB_PATH=/opt/local/lib/mpich-mp
```

Verify the inputs used by `configure`:

```bash
test -f "$RMPI_INCLUDE/mpi.h"
test -f "$RMPI_LIB_PATH/libmpi.dylib"
```

Do not set `RMPI_LIBS`; npRmpi's configure script does not use that variable.
Do not set `CC=mpicc` or `CXX=mpicxx` as part of the ordinary recipe. The MPI
include and library paths are supplied explicitly, while compilation should
remain aligned with R's configured toolchain.

Inspect that toolchain before building:

```bash
R CMD config CC
R CMD config CXX
R CMD config FC
R CMD config FLIBS
```

On Apple Silicon framework R 4.6, the validated `FLIBS` value is:

```text
/Library/Frameworks/R.framework/Resources/lib/libgfortran.5.dylib /Library/Frameworks/R.framework/Resources/lib/libquadmath.0.dylib
```

If `FLIBS` instead refers to missing or stale compiler paths, repair the local
R toolchain or user `~/.R/Makevars` before building. Do not hide missing
libraries with ad hoc dynamic-library path changes.

## 3. Build, install, and check the tarball

Run from the repository root:

```bash
VERSION=$(awk -F': *' '/^Version:/{print $2; exit}' DESCRIPTION)
R CMD build .
R CMD INSTALL "npRmpi_${VERSION}.tar.gz"
R CMD check --as-cran "npRmpi_${VERSION}.tar.gz"
```

For release work, install the tarball into a clean private library and test
that installed build. A successful source-loaded or in-tree run is not
installed-package proof.

On macOS, verify that the installed shared library resolves to MPICH and not
Open MPI:

```bash
NPRMPI_SO=$(Rscript -e 'cat(system.file("libs", paste0("npRmpi", .Platform$dynlib.ext), package = "npRmpi"))')
otool -L "$NPRMPI_SO"
```

The output should contain `/opt/local/lib/mpich-mp/libmpi` and must not contain
an Open MPI library path.

## 4. Installed runtime smoke

Run this in a fresh R process after installing from the tarball:

```r
smoke <- function() {
  library(npRmpi)
  options(np.messages = FALSE)

  npRmpi.init(nslaves = 1, quiet = TRUE)
  on.exit(try(npRmpi.quit(force = TRUE), silent = TRUE), add = TRUE)

  set.seed(1)
  x <- runif(200)
  y <- sin(2 * pi * x) + rnorm(200, sd = 0.2)
  bw <- npregbw(y ~ x, regtype = "ll", bwmethod = "cv.ls")
  fit <- npreg(bws = bw)
  stopifnot(length(fitted(fit)) == length(y), all(is.finite(fitted(fit))))

  cat("NPRMPI_MPICH_SMOKE_PASS\n")
}

smoke()
```

On macOS, the package defaults to retaining spawned workers for reuse.
`npRmpi.quit(force = TRUE)` is therefore required when the smoke is finished
and clean teardown is part of the evidence.

## 5. Batch/cluster attach mode

Launch a pre-created MPI world with:

```bash
mpiexec -n 4 Rscript foo.R
```

Use this structure in `foo.R`:

```r
library(npRmpi)
is_master <- npRmpi.init(mode = "attach", autodispatch = TRUE,
                         np.messages = FALSE)

if (isTRUE(is_master)) {
  # ... np* calls ...
  npRmpi.quit(mode = "attach")
}
```

Use a fresh launch for each independent attach-mode job. Confirm that no
route-scoped `mpiexec`, HYDRA, R, or slave processes remain afterward.

## 6. Network-interface diagnostics

MPICH's OFI transport may require an explicit interface on some hosts. Treat
interface selection as a host-specific diagnostic, not a package default.
First identify the intended interface, then apply the setting to one command:

```bash
FI_TCP_IFACE=en0 Rscript your-smoke.R
```

For a loopback-only diagnostic:

```bash
FI_TCP_IFACE=lo0 FI_PROVIDER=tcp FI_SOCKETS_IFACE=lo0 Rscript your-smoke.R
```

Do not use a loopback fallback as release acceptance for multi-host MPI.

`NP_RMPI_SKIP_INIT=1` deliberately prevents MPI initialization while loading
the package. It may help localize restricted check infrastructure, but it is
not runtime validation and must not be used for a positive release claim.

## 7. Release evidence

A release claim requires more than `R CMD check`. At minimum retain:

- exact source-tarball checksum and source head;
- configure, compile, install, and `otool -L` linkage proof;
- installed session and attach/profile tests at the required worker counts;
- numerical parity with serial np where mathematically expected;
- lifecycle, force-close, relaunch, and orphan-process proof;
- categorized examples, tests, demos, benchmarks, and rchk results;
- large-workload scaling evidence.

The known macOS MPICH teardown exit 137 is ignorable only when the proof log
contains the substantive pass marker/table and a subsequent process scan shows
that no route-scoped workers remain.
