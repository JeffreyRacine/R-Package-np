# npRmpi Demo Run Guide

This directory is the canonical demo harness for `npRmpi` timing/behavior
checks. It supports quick smoke runs and cross-version regression comparisons
against archived CRAN/candidate demo output.

## What This Runs

- `serial` mode: `np` scripts (`*_serial.R`) via `R CMD BATCH --no-save`
- `session` mode: `npRmpi` scripts (`*_npRmpi_session.R`) that call
  `npRmpi.init(nslaves=...)` from one master R process
- `attach` mode: `npRmpi` attach scripts (`*_npRmpi_attach.R`) via `mpiexec ... Rscript --no-save`; startup entry point is `npRmpi.init(mode="attach")`
- `profile` mode: `npRmpi` manual-broadcast scripts (`*_npRmpi_profile.R`) via `mpiexec ... R CMD BATCH --no-save`

`runall` executes in this order:
1. `serial`
2. `session` for each value in `SESSION_SLAVES`
3. for each value in `MPI_RANKS`:
   - `attach`
   - `profile`
4. parse transcripts into `timing/demo_results.csv` and text/TeX summaries

Default compute tiers:

- `SESSION_SLAVES="1 2 3"`
- `MPI_RANKS="2 3 4"`

Remember that `MPI_RANKS=4` means master plus three workers, so it is comparable
to `SESSION_SLAVES=3`.

## Prerequisites

1. Install `np` and `npRmpi` first (matching versions for your run).
2. MPI launcher available (`mpiexec`/`mpirun`).
3. Keep the expected folder layout (`demo/serial`, `demo/n_2_attach`, `demo/n_2_profile`, etc.).

## Core Files

- `runall`: orchestrates full matrix
- `makefile`: mode-specific launcher logic
- `../inst/demo_utils.R`: sample-size and machine-readable result helpers
- `parse_demo_results.R`: transcript parser used by `timing`
- `*.R`: demo scripts
- `../inst/Rprofile`: canonical profile startup (also available as `system.file("Rprofile", package="npRmpi")`)

The `makefile` is the source of truth for launch semantics:
- `attach`: timeout + cleared profile envs (`R_PROFILE_USER`, `R_PROFILE`) + optional `FI_*` env overrides
- `profile`: timeout + explicit `R_PROFILE_USER` + cleared `R_PROFILE` + optional `FI_*` env overrides + `NP_RMPI_PROFILE_RECV_TIMEOUT_SEC`
- all mode loops (`serial`, `attach`, `profile`) are fail-fast per demo; any failed demo exits non-zero immediately (no masked failures).
- session scripts run under `setsid` when available so MPI spawn/finalize
  signals do not terminate the launcher shell; the launcher cleans up spawned
  slave daemons after each session demo.
- attach demo scripts execute estimator bodies on master rank only (`mpi.comm.rank(0L) == 0L`) and finalize with `npRmpi.quit(mode="attach", ...)`.
- profile demo scripts do not call `npRmpi.init()`; startup comes from `inst/Rprofile`, then scripts use `np.mpi.initialize()` plus explicit `mpi.bcast.*` calls.

Profile startup contract (required):
- provide exactly one profile source per profile run:
  - either local `.Rprofile` in working directory, or
  - `R_PROFILE_USER=<profile-path>` (used by `makefile`);
- do not export both `R_PROFILE_USER` and `R_PROFILE` to the same file in one launch;
- `R CMD BATCH --no-save` is supported for profile mode when this contract is respected.
- package-level guard (`inst/Rprofile`) now hard-fails on dual-source profile startup with a remediation message.
- launchers export `NP_RMPI_PROFILE_RECV_TIMEOUT_SEC=$(TIMEOUT_SEC)` so blocked profile receives fail-fast.

## Tiny Fast Smoke

Use `n=100`, `NP=2` (`nslaves=1`) for quick validation.

### Serial

```bash
mkdir -p /Users/jracine/Development/np-npRmpi/demo/serial
cd /Users/jracine/Development/np-npRmpi/demo/serial
NP_DEMO_N=100 make -f ../makefile MODE=serial
```

### Session

```bash
mkdir -p /Users/jracine/Development/np-npRmpi/demo/session_direct
cd /Users/jracine/Development/np-npRmpi/demo/session_direct
NP_DEMO_N=100 make -f ../makefile MODE=session NSLAVES=1 DEMOS=npcdensls
```

### Attach

```bash
mkdir -p /Users/jracine/Development/np-npRmpi/demo/n_2_attach
cd /Users/jracine/Development/np-npRmpi/demo/n_2_attach
NP_DEMO_N=100 make -f ../makefile MODE=attach NP=2
```

### Profile

```bash
mkdir -p /Users/jracine/Development/np-npRmpi/demo/n_2_profile
cd /Users/jracine/Development/np-npRmpi/demo/n_2_profile
NP_DEMO_N=100 make -f ../makefile MODE=profile NP=2
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

Fractional calibrated-size smoke:

```bash
cd /Users/jracine/Development/np-npRmpi/demo
NP_DEMO_N_FRAC=0.1 ./runall
```

One-family four-route pilot:

```bash
cd /Users/jracine/Development/np-npRmpi/demo
DEMO_SET=pilot NP_DEMO_N=100 SESSION_SLAVES=1 MPI_RANKS=2 ./runall
```

Demo sets:

- `DEMO_SET=all`: makefile default demo list
- `DEMO_SET=pilot`: narrow `npcdensls` four-route pilot
- `DEMO_SET=conditional-density`: `npcdensls` and `npcdensml`
- `DEMO_SET=conditional-core`: `npcdensls`, `npcdensml`, and `npcdistls`
- `DEMO_SET=core-scaling`: initial workhorse subset using existing scripts
- `DEMO_SET=nomad`: reserved until explicit NOMAD demo scripts are added

You can always override the set with `DEMOS="npcdensls npreglcls"`.

Cross-version archive run:

```bash
cd /Users/jracine/Development/np-npRmpi/demo
RESULTS_ROOT=/Users/jracine/Development/demo_archives \
RUN_ID=npRmpi_0.70-2_candidate_$(date +%Y%m%d_%H%M%S) \
NP_DEMO_N=100 ./runall
```

Reset demo folder to baseline state after a run:

```bash
cd /Users/jracine/Development/np-npRmpi/demo
./cleanup
```

Shortcut:

```bash
cd /Users/jracine/Development/np-npRmpi/demo
./runall --cleanup
```

## Output Locations

- Default run outputs: `results/<run-id>/...`
- Serial outputs: `results/<run-id>/serial/*.Rout`
- Session outputs: `results/<run-id>/session/slaves_01/*.Rout`, etc.
- Attach outputs: `results/<run-id>/mpi_launch/ranks_02/attach/*.Rout`, etc.
- Profile outputs: `results/<run-id>/mpi_launch/ranks_02/profile/*.Rout`, etc.
- Parsed output: `results/<run-id>/timing/demo_results.csv`
- Cleanup helper: `./cleanup` removes run-generated folders/files and restores tracked timing outputs when run inside the repo

The parser also tolerates legacy direct folders such as `serial/`,
`n_2_attach/`, and `n_2_profile/`.

## Execution Context Rules (Must Follow)

### A) Running inside the repo

Use these paths:
- root: `/Users/jracine/Development/np-npRmpi/demo`
- attach/profile working dirs: `n_2_attach`, `n_2_profile`, etc.
- canonical profile path: `../inst/Rprofile` (resolved by `makefile`)

Recommended:

```bash
cd /Users/jracine/Development/np-npRmpi/demo
./runall
```

### B) Running from a copied demo folder outside the repo

You must preserve relative layout:
- copied root contains `makefile`, `runall`, and `*.R`
- run from subdirs (`serial`, `n_2_attach`, `n_2_profile`) with `make -f ../makefile ...`

For profile mode, set `RPROFILE` explicitly (do not rely on accidental relative matches):

```bash
RPROFILE=$(Rscript --no-save -e 'cat(system.file("Rprofile", package="npRmpi"))')
cd /path/to/copied/demo/n_2_profile
make -f ../makefile MODE=profile NP=2 NP_DEMO_N=100 RPROFILE="$RPROFILE"
```

If you also copy a local profile file, point `RPROFILE` to that explicit absolute path.

Optional MPI/libfabric network overrides (only when needed on your host):

```bash
FI_PROVIDER=tcp FI_TCP_IFACE=en0 FI_SOCKETS_IFACE=en0 \
make -f ../makefile MODE=attach NP=2 NP_DEMO_N=100
```

By default, the launcher does not force `FI_*`; it uses host MPI defaults.

## Troubleshooting

### Profile error: `could not find function "mpi.bcast.cmd"`

This indicates profile startup was not applied to ranks.
1. Confirm effective command:

```bash
cd /path/to/demo/n_2_profile
make -f ../makefile MODE=profile NP=2 NP_DEMO_N=100 DEMOS='npcdensls' -n run-profile
```

2. Verify the printed command contains `-env R_PROFILE_USER <expected path>`.
3. Re-run with explicit `RPROFILE=...` absolute path.
4. Ensure shell/session does not also export `R_PROFILE` to that same path.
5. Optional: set `NP_RMPI_PROFILE_RECV_TIMEOUT_SEC=<seconds>` to fail-fast if workers block waiting for broadcast commands.

### Attach/profile appears hung

1. Test with tiny smoke first (`NP_DEMO_N=100`, `NP=2`).
2. Run one script at a time (`DEMOS=<name>` in `make`).
3. Monitor workers directly:

```bash
pgrep -f 'npRmpi_attach.R|npRmpi_profile.R' | paste -sd, -
ps -p '<pid1,pid2>' -o 'pid,stat,etime,%cpu,command'
```

4. If one worker stays near `0%` while another stays high for an extended period, abort and retry.
5. Clean stale daemons before rerun:

```bash
pkill -f slavedaemon.R || true
pkill -f Rslaves.sh || true
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
