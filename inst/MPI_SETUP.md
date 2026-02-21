# MPI Setup (Current Workflow)

`npRmpi` supports two modern startup modes via `npRmpi.init()`.

## 1) Interactive (spawn mode)

Use when working inside an interactive R session and you want `npRmpi` to spawn
slaves for you.

```r
library(npRmpi)
npRmpi.init(nslaves = 1, mode = "spawn", autodispatch = TRUE, np.messages = FALSE)

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
npRmpi.init(mode = "attach", autodispatch = TRUE, np.messages = FALSE)

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
