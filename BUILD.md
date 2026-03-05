# Build (np)

## Build Tarball

```bash
cd /Users/jracine/Development
R CMD build np-master
```

This produces `np_0.70-0.tar.gz` in `/Users/jracine/Development`.

## Install

```bash
cd /Users/jracine/Development
R CMD INSTALL np_0.70-0.tar.gz
```

## Quick Load Check

```bash
R -q -e 'library(np); sessionInfo()'
```

## Check

```bash
cd /Users/jracine/Development
R CMD check --as-cran np_0.70-0.tar.gz
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
