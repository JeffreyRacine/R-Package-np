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
4. parse transcripts into `timing/demo_results.csv`, `demo_results_wide.csv`,
   `timing_all.dat`, and `demo_summary.qmd`
5. render `demo_summary.html` when Quarto is available

Default compute tiers:

- `SESSION_SLAVES="1 2 3"`
- `MPI_RANKS="2 3 4"`

Remember that `MPI_RANKS=4` means master plus three workers, so it is comparable
to `SESSION_SLAVES=3`.

## Prerequisites

1. Install `np` and `npRmpi` first (matching versions for your run).
2. MPI launcher available (`mpiexec`/`mpirun`).
3. Keep the expected folder layout when running copied demo folders; the
   canonical in-repo runner writes under `demo/results/<run-id>/`.

## Core Files

- `runall`: top-level compatibility shim for `tools/runall`
- `makefile`: top-level compatibility include for `tools/makefile`
- `tools/runall`: orchestrates full matrix
- `tools/makefile`: mode-specific launcher logic
- `tools/monitor_demo_cpu.sh`: optional process-level CPU sidecar for MPI runs
- `tools/summarize_cpu_monitor.R`: summarizes the CPU sidecar and flags low-utilization stages
- `tools/timing`: parser/report wrapper
- `tools/parse_demo_results.R`: transcript parser and Quarto report generator
- `../inst/demo_utils.R`: sample-size and machine-readable result helpers
- `../inst/demo_matrices/*.csv`: family/case matrices for modern demos
- `../inst/demo_family_*.R`: shared family runners used by thin mode scripts
- `*.R`: demo scripts
- `../inst/Rprofile`: canonical profile startup (also available as `system.file("Rprofile", package="npRmpi")`)

The `tools/makefile` is the source of truth for launch semantics:
- `attach`: timeout + cleared profile envs (`R_PROFILE_USER`, `R_PROFILE`) + optional `FI_*` env overrides
- `profile`: timeout + explicit `R_PROFILE_USER` + cleared `R_PROFILE` + optional `FI_*` env overrides + `NP_RMPI_PROFILE_RECV_TIMEOUT_SEC`
- all mode loops are fail-fast per demo by default; set
  `CONTINUE_ON_ERROR=true` for archive/sentinel comparison runs where known
  failures should be recorded and the remaining demos should continue.
- session scripts run under `setsid` when available so MPI spawn/finalize
  signals do not terminate the launcher shell; the launcher cleans up spawned
  slave daemons after each session demo.
- session mode skips demos that do not yet have an explicit
  `*_npRmpi_session.R` script, so broad legacy demo sets remain runnable while
  session coverage is expanded family by family.
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
- MPI launch timeout defaults to `TIMEOUT_SEC=600` seconds. Increase this
  explicitly for deliberately long diagnostics; do not use multi-hour timeouts
  for routine sentinel runs.
- set `CPU_MONITOR=true` to write `cpu_monitor.csv` and
  `cpu_monitor_summary.csv` under the run directory. Stages with rank median
  CPU below roughly `80-90%` should be treated as suspicious enough to inspect
  when the corresponding timing row is long enough to be meaningful.

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
NP_DEMO_N=100 make -f ../makefile MODE=session NSLAVES=1 DEMOS=npcdens
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

Full matrix intended for archive comparison:

```bash
cd /Users/jracine/Development/np-npRmpi/demo
NP_DEMO_TIER=sentinel CONTINUE_ON_ERROR=true CPU_MONITOR=true TIMEOUT_SEC=600 RUNALL_REQUIRE_FULL_MATRIX=true ./runall
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
- `DEMO_SET=pilot`: narrow matrix-driven `npcdens` four-route pilot
- `DEMO_SET=conditional-density`: matrix-driven `npcdens`
- `DEMO_SET=conditional-core`: matrix-driven `npcdens` and `npcdist`
- `DEMO_SET=core-scaling`: converted core subset currently covering
  matrix-driven `npreg`, matrix-driven `npcdens`, matrix-driven `npcdist`,
  matrix-driven `npudens`, and matrix-driven `npudist`
- `DEMO_SET=nomad`: matrix-driven `npreg`, `npcdens`, `npcdist`, `npindex`, `npscoef`, and `npplreg` NOMAD smoke rows
- `DEMO_SET=semiparametric`: matrix-driven `npindex`, `npscoef`, and `npplreg`
- `DEMO_SET=unconditional`: matrix-driven `npudens` and `npudist`
- `DEMO_SET=tests`: matrix-driven statistical test demos:
  `npdeneqtest`, `npdeptest`, `npsdeptest`, `npsigtest`, `npsymtest`,
  and `npunitest`
- `DEMO_SET=auxiliary`: matrix-driven auxiliary demos:
  `npcmstest`, `npconmode`, `npcopula`, `npqreg`, and `npregiv`

You can always override the set with `DEMOS="npcdens npreg"`.

Matrix-driven family demos, beginning with `npreg`, accept:

- `NP_DEMO_TIER=smoke` or `NP_DEMO_TIER=sentinel`
- `NP_DEMO_CASES="npreg_lc_ls npreg_lp_nomad"` for selected rows
- `NP_DEMO_MATRIX=/absolute/path/to/custom.csv` for a custom matrix

The `npreg` smoke matrix currently covers `lc`, `ll`, explicit `lp` degree 1,
`cv.aic`, and NOMAD local-polynomial degree search. The heavier sentinel
matrix uses larger calibrated defaults for cross-version regression checks.
Sentinel `default_n` values are calibrated from current timings: rows that
already track roughly `serial / mpi_size` are left alone even when short, while
functioning but underpowered rows are increased. Suspicious flat or failing
rows are not hidden by increasing `n`; keep them visible for debugging.

The `npcdens` smoke matrix currently covers `lc` least-squares, `lc`
maximum-likelihood, `ll` least-squares, explicit `lp` degree 1, and NOMAD
local-polynomial degree search. The sentinel matrix preserves the previous
`npcdensls`/`npcdensml` calibrated defaults while adding the newer `ll`, `lp`,
and NOMAD surfaces.

The `npcdist` smoke matrix currently covers `lc` least-squares, `ll`
least-squares, explicit `lp` degree 1, and NOMAD local-polynomial degree
search. The sentinel matrix preserves the previous `npcdistls` calibrated
default while adding the newer `ll`, `lp`, and NOMAD surfaces.

The `npindex` smoke matrix replaces the older Ichimura/Klein-Spady split
scripts with rows for Ichimura `lc`, `ll`, explicit `lp` degree 1, Ichimura
LP NOMAD degree search, and the Klein-Spady binary route. It records beta and
bandwidth fields as regression-litmus sentinels.

The `npscoef` smoke matrix preserves the original smooth-coefficient data
generator while adding `ll`, explicit `lp` degree 1, and LP NOMAD degree-search
rows. It records objective, bandwidth, fitted-length, and evaluation-count
fields so the semiparametric route can be compared across serial, session,
attach, and profile/manual-broadcast modes.

The `npplreg` smoke matrix keeps the original partially linear data generator
and adds local-linear, explicit degree-1 local-polynomial, and LP NOMAD
degree-search rows. Because `npplreg` stores multiple child bandwidths, the
timing parser records the first two flattened bandwidth and coefficient
sentinels for regression-litmus comparisons.

The `npudens` smoke matrix replaces the older `npudensls`/`npudensml` split
with one unconditional density demo covering both least-squares and likelihood
cross-validation. Each row now runs both bandwidth selection and `npudens()`,
recording objective, bandwidth, fitted length, and the first density estimate.

The `npudist` smoke matrix replaces `npudistcdf` with one unconditional
distribution demo covering CDF cross-validation. It runs both `npudistbw()` and
`npudist()`, recording objective, bandwidth, fitted length, and the first CDF
estimate.

The statistical test demos now share a matrix-driven runner with calibrated
sample-size and bootstrap-count rows. Smoke matrices use the minimum supported
bootstrap count, while sentinel matrices preserve the older heavier defaults.
The parser records representative numeric test statistics and p-values when
available.

The auxiliary demos now share a matrix-driven runner with deterministic
simulated data. This removes fixed-size real-data dependencies from the scaling
surface for `npcmstest`, `npconmode`, and `npqreg`, while preserving the public
workflow each demo is meant to exercise.

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
- Wide parsed output: `results/<run-id>/timing/demo_results_wide.csv`
- Human-readable terminal summary: `results/<run-id>/timing/timing_all.dat`
- Planned run inventory: `results/<run-id>/RUN_INVENTORY.tsv`
- Observed mode/rank/slave ledger: `results/<run-id>/RUN_OBSERVED_MATRIX.tsv`
- Matrix validation status: `results/<run-id>/RUN_MATRIX_VALIDATION.tsv`
- Top-level run status: `results/<run-id>/RUN_STATUS.tsv`
- Quarto report source: `results/<run-id>/timing/demo_summary.qmd`
- Optional rendered report, when Quarto is installed:
  `results/<run-id>/timing/demo_summary.html`
- Cleanup helper: `./cleanup` removes run-generated folders/files when run
  inside the repo

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
- copied root contains `makefile`, `runall`, `tools/`, and `*.R`
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
make -f ../makefile MODE=profile NP=2 NP_DEMO_N=100 DEMOS='npcdens' -n run-profile
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
