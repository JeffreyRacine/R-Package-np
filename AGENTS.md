# AGENTS.md

Use `/Users/jracine/Development/AGENTS.md` as the canonical workflow and benchmarking policy.
Hard default: enforce `Sweep Safety Protocol (Default Across Scope Repos)` from the canonical AGENTS file.

Repo-specific note:
- Merge (do not overwrite) R-layer and man-page files when porting from `np-master`.
- Any new or ported benchmark/function must follow MPI protocol from `/Users/jracine/Development/AGENTS.md`: start/stop slaves, broadcast commands/options/data explicitly, support passable `nslaves`/`rslaves`, and keep parity with serial seed/DGP/options.
- Performance evidence policy for this repo follows canonical AGENTS exactly:
  - `times=5` is smoke only; decision screening starts at interleaved paired-seed `times=25`;
  - predeclare MEI (absolute seconds + relative percent), use paired deltas + CI, and allow early accept at `times=25` only for obvious wins;
  - otherwise confirm with interleaved paired-seed `times>=100` (or higher when pilot variance implies it);
  - accept on performance grounds only when paired mean and paired median agree;
  - all gated methods/modes must agree in direction, otherwise anomaly rationalization versus MEI is mandatory;
  - default MPI performance gates to `nslaves=1`; run `nslaves=0` only for concrete master-only diagnostics.
- Keep shared regression summary telemetry in sync with `np-master`: `num.feval`, `num.feval.fast`, `num.feval.fallback`.
- Keep `man/np.kernels.Rd` and kernel/options/plot cross-links aligned, while preserving package-specific `npRmpi` docs text.
- Runtime contract: `npRmpi` must never call serial `np` estimator/plot/bootstrap/helper code paths and must not silently fall back to serial execution when MPI mode is selected.
- Demo execution contract (`np-npRmpi/demo`):
  - treat `demo/makefile` as launcher source of truth for attach/profile;
  - run from `demo` subdirectories (`serial`, `n_*_attach`, `n_*_profile`) using `make -f ../makefile ...`;
  - in-repo profile path should resolve to `../inst/Rprofile`; for copied demo folders outside repo, set `RPROFILE` explicitly to an absolute profile path (typically `system.file("Rprofile", package="npRmpi")`).
  - for profile/manual-broadcast launch, provide exactly one profile source (`.Rprofile` in cwd or `R_PROFILE_USER`); never set both `R_PROFILE_USER` and `R_PROFILE` to the same file in a single `mpiexec` launch.
  - `R CMD BATCH --no-save` is valid for profile mode when the single-source profile rule is followed.
  - profile startup guard in `inst/Rprofile` is expected to fail-fast on dual-source profile wiring with a remediation message.
- Demo failure recovery contract:
  - reproduce first with `NP_DEMO_N=100`, `NP=2`, and `DEMOS='<single demo>'`;
  - validate effective profile command with `make -n run-profile` and confirm `R_PROFILE_USER` path;
  - for suspected hangs, monitor worker split with `ps -p '<pid1,pid2>' -o 'pid,stat,etime,%cpu,command'`;
  - if one worker remains near idle while another remains busy, abort and clean stale `slavedaemon.R`/`Rslaves.sh` before retry.
