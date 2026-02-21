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
