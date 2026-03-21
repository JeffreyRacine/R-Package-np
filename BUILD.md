# Build (np)

## One-Command Preflight

```bash
cd /Users/jracine/Development
./release_preflight.sh np
```

This builds `np`, runs the shared package/gallery sync audit, and then runs the
final tarball-first `R CMD check --as-cran`.

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

## Check

```bash
cd /Users/jracine/Development
./package_gallery_sync_audit.sh
VERSION=$(awk -F': *' '/^Version:/{print $2; exit}' np-master/DESCRIPTION)
R CMD check --as-cran "np_${VERSION}.tar.gz"
```

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
