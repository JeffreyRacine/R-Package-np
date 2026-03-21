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

## One-Command Preflight

```bash
cd /Users/jracine/Development
./release_preflight.sh npRmpi
```

This builds `npRmpi`, runs the shared package/gallery sync audit, and then runs
the final tarball-first `R CMD check --as-cran`.

## Build Tarball

```bash
cd /Users/jracine/Development
R CMD build np-npRmpi
VERSION=$(awk -F': *' '/^Version:/{print $2; exit}' np-npRmpi/DESCRIPTION)
```

This produces `npRmpi_${VERSION}.tar.gz` in `/Users/jracine/Development`.

## Install

```bash
cd /Users/jracine/Development
VERSION=$(awk -F': *' '/^Version:/{print $2; exit}' np-npRmpi/DESCRIPTION)
R CMD INSTALL "npRmpi_${VERSION}.tar.gz"
```

## Quick Load Check

```bash
R -q -e 'library(npRmpi); sessionInfo()'
```

## Check

```bash
cd /Users/jracine/Development
./package_gallery_sync_audit.sh
VERSION=$(awk -F': *' '/^Version:/{print $2; exit}' np-npRmpi/DESCRIPTION)
R CMD check --as-cran "npRmpi_${VERSION}.tar.gz"
```

## Release-Surface Audit

When vignette names, startup routing, package help routing, gallery package
links, or source-tarball vignette packaging change, run the shared audit before
signoff:

```bash
cd /Users/jracine/Development
./release_preflight.sh npRmpi
```

For `npRmpi`, this audit also checks that the source tarball keeps
`build/vignette.rds`.

If gallery routing pages also changed, use:

```bash
cd /Users/jracine/Development
./release_preflight.sh --render-gallery npRmpi
```

### MPI Example Modes During Check

- Default `R CMD check --as-cran` behavior:
  - MPI spawn examples are skipped in check context to avoid MPICH teardown
    killing the parent check process on some systems.
- To force MPI spawn examples to run during check, use the exact sequence below.

### Forced MPI Example Check (Exact Sequence)

```bash
cd /Users/jracine/Development/np-npRmpi/man
./run

cd /Users/jracine/Development
R CMD build np-npRmpi

cd /Users/jracine/Development
VERSION=$(awk -F': *' '/^Version:/{print $2; exit}' np-npRmpi/DESCRIPTION)
NP_RMPI_RUN_MPI_EXAMPLES_IN_CHECK=1 R CMD check --as-cran "npRmpi_${VERSION}.tar.gz"

cd /Users/jracine/Development/np-npRmpi/man
./dontrun
```

## Runtime Modes

- Interactive R session:
  - `npRmpi.init(nslaves=...)` for session mode (the `spawn` code path)
- Cluster/batch under `mpiexec`:
  - start script with `npRmpi.init(mode="attach", autodispatch=TRUE, np.messages=FALSE)`
  - end script with `npRmpi.quit(mode="attach")`
- Manual-broadcast/profile mode:
  - start ranks with `inst/Rprofile` (or `R_PROFILE_USER`) and explicit `mpi.bcast.*` calls
  - use `np.mpi.initialize()` rather than `npRmpi.init()` as the workflow initializer

## Canonical Implementation Directive (2026-03-05)

This repository follows a strict canonical execution rule:

1. One canonical implementation per method (outside explicit `np.tree` branching).
2. Unsupported configurations must fail fast with explicit `stop(...)` diagnostics.
3. No silent remap/coercion of user-selected options (for example `bwmethod`, `regtype`, kernels, `cv.iterate`, or bounds transforms).
4. No hidden alternate execution paths for the same method semantics.
5. All fit-defining options (for example `degree`, `basis`, `bernstein.basis`, kernels, and bounds) must be propagated and used by the canonical path.
6. `np.tree=FALSE` is the default; when `np.tree=TRUE`, behavior must remain semantics-preserving and option-compatible with the canonical path.
7. Remove or reject legacy/debug compatibility branches that add redundant runtime overhead once canonical behavior exists.
