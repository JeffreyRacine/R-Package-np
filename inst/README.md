# Installed Runtime Assets (`inst/`)

This directory should contain only files that are installed with the package
and needed at runtime (or package metadata such as citation/legal text).

## Runtime scripts

- `slavedaemon.R`: worker-loop script used by `mpi.spawn.Rslaves()`.
- `Rslaves.sh`: POSIX launcher used by `mpi.spawn.Rslaves()`.
- `MacR64slaves.sh`: legacy Intel-mac launcher retained for compatibility.

## Metadata

- `CITATION`: citation metadata returned by `citation("npRmpi")`.
- `COPYRIGHTS`: copyright notice.

## User guidance

- `MPI_SETUP.md`: current, supported startup patterns for interactive and
  cluster/batch use.

Historical build/install notes and one-off development artifacts were removed
from `inst/` to keep installed payload minimal and current.

## Canonical Implementation Directive (2026-03-05)

This repository follows a strict canonical execution rule:

1. One canonical implementation per method (outside explicit `np.tree` branching).
2. Unsupported configurations must fail fast with explicit `stop(...)` diagnostics.
3. No silent remap/coercion of user-selected options (for example `bwmethod`, `regtype`, kernels, `cv.iterate`, or bounds transforms).
4. No hidden alternate execution paths for the same method semantics.
5. All fit-defining options (for example `degree`, `basis`, `bernstein.basis`, kernels, and bounds) must be propagated and used by the canonical path.
6. `np.tree=FALSE` is the default; when `np.tree=TRUE`, behavior must remain semantics-preserving and option-compatible with the canonical path.
7. Remove or reject legacy/debug compatibility branches that add redundant runtime overhead once canonical behavior exists.
