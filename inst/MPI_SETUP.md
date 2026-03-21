# MPI Setup (Current Workflow)

`npRmpi` supports two modern startup modes via `npRmpi.init()`.

## 1) Interactive (spawn mode)

Use when working inside an interactive R session and you want `npRmpi` to spawn
slaves for you.

```r
library(npRmpi)
npRmpi.init(nslaves = 1)

# ... run np* calls ...

npRmpi.quit()
```

## 2) Cluster/Batch (attach mode)

Use when MPI world is pre-launched (e.g. with `mpiexec`).

```bash
mpiexec -n 128 Rscript foo.R
```

Inside `foo.R`:

```r
library(npRmpi)
npRmpi.init(mode = "attach")

# ... run np* calls ...

npRmpi.quit(mode = "attach")
```

## Notes

- `.Rprofile` bootstrap files are not required for the supported
  `npRmpi.init()` workflow.
- `mode = "auto"` selects `attach` when `mpi.comm.size(0) > 1`, otherwise
  `spawn`.
- On some systems, explicitly setting MPI interface can help, e.g.
  `FI_TCP_IFACE=en0` (fallback: `FI_TCP_IFACE=lo0`).

## Canonical Implementation Directive (2026-03-05)

This repository follows a strict canonical execution rule:

1. One canonical implementation per method (outside explicit `np.tree` branching).
2. Unsupported configurations must fail fast with explicit `stop(...)` diagnostics.
3. No silent remap/coercion of user-selected options (for example `bwmethod`, `regtype`, kernels, `cv.iterate`, or bounds transforms).
4. No hidden alternate execution paths for the same method semantics.
5. All fit-defining options (for example `degree`, `basis`, `bernstein.basis`, kernels, and bounds) must be propagated and used by the canonical path.
6. `np.tree=FALSE` is the default; when `np.tree=TRUE`, behavior must remain semantics-preserving and option-compatible with the canonical path.
7. Remove or reject legacy/debug compatibility branches that add redundant runtime overhead once canonical behavior exists.
