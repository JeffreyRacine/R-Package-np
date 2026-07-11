# Building and validating npRmpi

The supported 0.70-5 release backend is MPICH. For macOS with MacPorts, follow
[BUILD_MPICH.md](BUILD_MPICH.md). Open MPI runtime validation is explicitly
waived for this release.

## Portable source-tarball workflow

Run from the repository root after configuring a supported MPI implementation:

```bash
VERSION=$(awk -F': *' '/^Version:/{print $2; exit}' DESCRIPTION)
R CMD build .
R CMD INSTALL "npRmpi_${VERSION}.tar.gz"
R CMD check --as-cran "npRmpi_${VERSION}.tar.gz"
```

Build from a clean repository. For release validation, install the generated
tarball into a clean private library and run tests against that installed
package rather than the live source tree.

## R compiler and LAPACK/BLAS linkage

MPI discovery is host-specific, but npRmpi links LAPACK, BLAS, and the Fortran
runtime through R's portable make variables:

```make
PKG_LIBS = @PKG_LIBS@ $(ARCHLIB) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)
```

Inspect the active R toolchain before building:

```bash
R CMD config CC
R CMD config CXX
R CMD config FC
R CMD config FLIBS
```

Do not replace R's compilers with MPI wrapper compilers merely to locate MPI.
Supply MPI type, include, and library paths through the supported configure
arguments or `RMPI_TYPE`, `RMPI_INCLUDE`, and `RMPI_LIB_PATH`.

## Installed linkage and startup proof

After installation, record the package and MPI runtime identities:

```bash
R -q -e 'library(npRmpi); print(npRmpi.session.info(comm = 0L)); sessionInfo()'
```

On systems with `otool` or `ldd`, verify that the installed `npRmpi` shared
library resolves to exactly the intended MPI implementation. A mixed-backend
build is a hard failure.

## Runtime modes

- Interactive/session mode: `npRmpi.init(nslaves = ...)` dynamically spawns a
  worker pool. On macOS, use `npRmpi.quit(force = TRUE)` when clean teardown is
  required; the non-force default may retain workers for reuse.
- Attach mode: prelaunch a world with `mpiexec`, call
  `npRmpi.init(mode = "attach")`, and close it with
  `npRmpi.quit(mode = "attach")` on the master.
- Profile/manual-broadcast mode: launch ranks with the installed `Rprofile`
  and use explicit `mpi.bcast.*` calls plus `np.mpi.initialize()`.

See the installed `MPI_SETUP.md` for examples.

## MPI examples during package checks

Ordinary `R CMD check --as-cran` skips spawn examples that can destabilize the
parent check process on some MPI/macOS combinations. To force those examples
for a dedicated MPI check, use the repository helpers and restore the normal
Rd state afterward:

```bash
(cd man && ./run)
VERSION=$(awk -F': *' '/^Version:/{print $2; exit}' DESCRIPTION)
R CMD build .
NP_RMPI_RUN_MPI_EXAMPLES_IN_CHECK=1 \
  R CMD check --as-cran "npRmpi_${VERSION}.tar.gz"
(cd man && ./dontrun)
```

Do not leave the live documentation tree in the forced-example state.

## Release-validation expectations

A CRAN-ready claim requires the enhanced installed-tarball protocol, including:

- clean source/head and tarball provenance;
- `R CMD check --as-cran`, examples, vignettes, and categorized tests;
- demos and configured benchmark sentinels;
- public estimator-family and S3-surface coverage;
- session/attach/profile lifecycle, cleanup, and relaunch proof;
- numerical equivalence with serial np where expected;
- large-workload MPI scaling at the required worker counts;
- rchk for native-code protection;
- reverse-dependency review when configured.

If the generic rchk container lacks MPI headers or build tooling, record an
infrastructure skip; do not describe that as native-code validation.

Timing reports must record `R.version.string`, the active BLAS, MPI version,
worker count, and backend linkage. A positive release claim must not rely on a
source-loaded probe or a diagnostic fallback.
