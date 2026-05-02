# Build (np)

## One-Command Preflight

```bash
cd /Users/jracine/Development
./release_preflight.sh np
```

This builds `np`, runs the shared package/gallery sync audit, and then runs the
final tarball-first `R CMD check --as-cran`.

For CRAN-ready release claims, use the hardened revdep-aware release gate
instead of relying on `release_preflight.sh` alone:

```bash
cd /Users/jracine/Development
./release_protocol/run_np_release_gate.sh
```

This gate builds from source, installs into a private library, runs installed
smokes, runs tarball-first `R CMD check --as-cran`, inventories CRAN reverse
dependencies, and can run the focused revdep checks required after the
2026-05-01 bandwidth-dispatch compatibility regression. It also runs the
containerized `rchk` native-code protection check when feasible
(`RUN_RCHK=auto` by default).

## LAPACK/BLAS Linkage

Unix builds must link through R's portable make variables, not host-specific
library paths:

```make
PKG_LIBS = $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)
```

On Apple Silicon framework R 4.6 hosts, verify the local R configuration before
release builds:

```bash
/Library/Frameworks/R.framework/Resources/bin/R CMD config FLIBS
```

For the current framework setup, `FLIBS` should resolve to:

```text
/Library/Frameworks/R.framework/Resources/lib/libgfortran.5.dylib /Library/Frameworks/R.framework/Resources/lib/libquadmath.0.dylib
```

If it instead points at stale `/opt/gfortran/...` paths, refresh the host
`Makeconf`/user `Makevars` as described in
`/Users/jracine/Development/SKILLS.md` before running release preflight.

## Build Tarball

```bash
cd /Users/jracine/Development
R CMD build np-master
VERSION=$(awk -F': *' '/^Version:/{print $2; exit}' np-master/DESCRIPTION)
```

This produces `np_${VERSION}.tar.gz` in `/Users/jracine/Development`.

## Install

```bash
cd /Users/jracine/Development
VERSION=$(awk -F': *' '/^Version:/{print $2; exit}' np-master/DESCRIPTION)
R CMD INSTALL "np_${VERSION}.tar.gz"
```

## Quick Load Check

```bash
R -q -e 'library(np); sessionInfo()'
```

## Timing Provenance

When collecting timing or regression artifacts, always record:

```bash
R -q -e 'cat("R.version.string =", R.version.string, "\n"); cat("BLAS =", sessionInfo()$BLAS, "\n")'
```

On macOS CRAN-binary R, `sessionInfo()$BLAS` reflects the active
`libRblas.dylib` choice. On Apple Silicon, that choice can materially affect
BLAS-heavy timings, so treat `R.version.string` and `sessionInfo()$BLAS` as
part of the benchmark environment.

## Check

```bash
cd /Users/jracine/Development
./package_gallery_sync_audit.sh
VERSION=$(awk -F': *' '/^Version:/{print $2; exit}' np-master/DESCRIPTION)
R CMD check --as-cran "np_${VERSION}.tar.gz"
```

If a change touches formula interfaces, `...` forwarding, bandwidth
constructors, compatibility dispatch helpers, exported defaults, warning/error
contracts, estimator semantics, or object structure, run the revdep-aware
release gate and focused revdep lane before submission.

If a change touches `src/`, registered native interfaces, or code that changes
the shape/lifetime of objects passed to `.C`, `.Call`, or `.Fortran`, require
local `rchk` proof when infrastructure is available:

```bash
cd /Users/jracine/Development
RUN_RCHK=1 ./release_protocol/run_np_release_gate.sh
```

Use `RUN_RCHK=auto` for ordinary full release rehearsal; it records a precise
`SKIP` when Docker/rchk infrastructure is unavailable.

## Release-Surface Audit

When vignette names, startup routing, package help routing, or gallery package
links change, run the shared audit before signoff:

```bash
cd /Users/jracine/Development
./release_preflight.sh np
```

If gallery routing pages also changed, use:

```bash
cd /Users/jracine/Development
./release_preflight.sh --render-gallery np
```

Then run the package check:

```bash
cd /Users/jracine/Development
./package_gallery_sync_audit.sh
VERSION=$(awk -F': *' '/^Version:/{print $2; exit}' np-master/DESCRIPTION)
R CMD check --as-cran "np_${VERSION}.tar.gz"
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
