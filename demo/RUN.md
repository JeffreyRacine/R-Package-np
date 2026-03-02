# npRmpi Demo Run Guide

This directory is the canonical demo harness for `npRmpi` timing/behavior checks.

## What This Runs

- `serial` mode: `np` scripts (`*_serial.R`) via `R CMD BATCH --no-save`
- `attach` mode: `npRmpi` attach scripts (`*_npRmpi_attach.R`) via `mpiexec ... Rscript --no-save`
- `profile` mode: `npRmpi` manual-broadcast scripts (`*_npRmpi_profile.R`) via `mpiexec ... R CMD BATCH --no-save`

`runall` executes in this order:
1. all `serial`
2. for each `NP in {2,3,4}`:
   - `attach`
   - `profile`

## Prerequisites

1. Install `np` and `npRmpi` first (matching versions for your run).
2. MPI launcher available (`mpiexec`/`mpirun`).
3. Run commands from this directory: `/Users/jracine/Development/np-npRmpi/demo`.

## Core Files

- `runall`: orchestrates full matrix
- `makefile`: mode-specific launcher logic
- `*.R`: demo scripts
- `../inst/Rprofile`: canonical profile startup (also available as `system.file("Rprofile", package="npRmpi")`)

Note:
- A local `.Rprofile` in this demo directory is optional.
- `makefile` profile mode passes `R_PROFILE_USER` through `mpiexec -env` so worker ranks receive the startup profile.

## Tiny Fast Smoke

Use `n=100`, `NP=2` (`nslaves=1`) for quick validation.

### Serial

```bash
cd /Users/jracine/Development/np-npRmpi/demo/serial
make -f ../makefile MODE=serial NP_DEMO_N=100
```

### Attach

```bash
cd /Users/jracine/Development/np-npRmpi/demo/n_2_attach
make -f ../makefile MODE=attach NP=2 NP_DEMO_N=100
```

### Profile

```bash
cd /Users/jracine/Development/np-npRmpi/demo/n_2_profile
make -f ../makefile MODE=profile NP=2 NP_DEMO_N=100
```

## Full Matrix

```bash
cd /Users/jracine/Development/np-npRmpi/demo
./runall
```

Optional smaller smoke:

```bash
cd /Users/jracine/Development/np-npRmpi/demo
NP_DEMO_N=100 ./runall
```

## Output Locations

- Serial outputs: `serial/*.Rout`
- Attach outputs: `n_2_attach/*.Rout`, `n_3_attach/*.Rout`, `n_4_attach/*.Rout`
- Profile outputs: `n_2_profile/*.Rout`, `n_3_profile/*.Rout`, `n_4_profile/*.Rout`

## Running From a Copied Demo Folder

If you copy this demo directory elsewhere, copy all of:
- `runall`
- `makefile`
- all `*.R` scripts

You do not need to copy `.Rprofile` if `npRmpi` is installed, since `makefile` can use:
- `Rscript -e 'cat(system.file(\"Rprofile\", package=\"npRmpi\"))'`

If needed, override explicitly:

```bash
RPROFILE=$(Rscript --no-save -e 'cat(system.file("Rprofile", package="npRmpi"))')
cd /path/to/copied/demo/n_2_profile
make -f ../makefile MODE=profile NP=2 NP_DEMO_N=100 RPROFILE=$RPROFILE
```

## Troubleshooting

### Profile error: `could not find function "mpi.bcast.cmd"`

This indicates profile startup was not applied to ranks.
Run profile mode through this `makefile` (which exports `R_PROFILE_USER` via `mpiexec -env`), or set `RPROFILE` explicitly as above.

### Attach/profile appears hung

1. Test with tiny smoke first (`NP_DEMO_N=100`, `NP=2`).
2. Run one script at a time (`DEMOS=<name>` in `make`).
3. Clean stale daemons before rerun:

```bash
pkill -f slavedaemon.R || true
pkill -f Rslaves.sh || true
```

