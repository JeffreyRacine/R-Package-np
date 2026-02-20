# npRmpi Callable Surface Inventory for Auto-Dispatch/Guardrails

## Scope
Inventory of user-callable R entry points in `/Users/jracine/Development/np-npRmpi` that must interoperate under an auto-dispatch helper (no manual `mpi.bcast.cmd(...)` required for normal usage), while avoiding hangs and clobbering.

Primary references:
- `/Users/jracine/Development/np-npRmpi/NAMESPACE`
- `/Users/jracine/Development/np-npRmpi/R/*.R`

## 1) Direct npRmpi User API (Top-Level Exports)
From `NAMESPACE` export block (lines 130-175), the core `np*` user-facing surface is:

1. Runtime/session:
- `np.mpi.initialize`
- `npRmpi.start`
- `npRmpi.stop`
- `npRmpi.session.info`

2. Estimation/bandwidth constructors:
- `npudensbw`, `npudistbw`, `npudens`, `npudist`
- `npregbw`, `npreg`
- `npplregbw`, `npplreg`
- `npcdensbw`, `npcdens`
- `npcdistbw`, `npcdist`
- `npindexbw`, `npindex`
- `npscoefbw`, `npscoef`
- `npqreg`, `npksum`, `npconmode`, `npcopula`
- `npregiv`, `npregivderiv`

3. Testing/utilities and other estimators:
- `npcmstest`, `npdeptest`, `npsdeptest`, `npdeneqtest`, `npqcmstest`, `npsigtest`, `npsymtest`, `npunitest`
- `npquantile`, `npseed`, `nptgauss`
- `npuniden.boundary`, `npuniden.reflect`, `npuniden.sc`
- `np.pairs`, `np.pairs.plot`, `uocquantile`

4. Shared helpers exposed to users:
- `npplot`
- `se`
- `gradients`

## 2) S3 Surface Registered in NAMESPACE
`NAMESPACE` registers a broad S3 surface (217 S3method entries). These are user-callable via generics and must be compatible with dispatch state.

Generic coverage counts:
- `print`: 26
- `summary`: 25
- `plot`: 17
- `predict`: 17
- `fitted`: 10
- `se`: 8
- `npplot`: 8
- `gradients`: 5
- `compute.bootstrap.errors`: 5
- `npqreg`, `npsigtest`: 5 each
- plus class constructors (`npreg*`, `npcdens*`, `npcdist*`, etc.) with 4-5 methods each.

## 3) Critical Transitive Paths (Must Be Included)
These are the highest-risk pathways for hangs when users do not manually broadcast:

1. `plot.<class>()` methods that delegate to `npplot(...)`:
- e.g. `plot.npregression`, `plot.condensity`, `plot.condistribution`, `plot.plregression`, `plot.qregression`, `plot.singleindex`, `plot.smoothcoefficient`, and bandwidth plot methods.

2. `predict.<class>()` methods that call back into `np*` constructors via `eval(...)`:
- `predict.npregression` -> `eval(npreg(...), envir=parent.frame())`
- `predict.condensity` -> `eval(npcdens(...), envir=parent.frame())`
- `predict.condistribution` -> `eval(npcdist(...), envir=parent.frame())`
- `predict.npdensity` -> `eval(npudens(...), envir=parent.frame())`
- `predict.npdistribution` -> `eval(npudist(...), envir=parent.frame())`
- `predict.plregression` -> `eval(npplreg(...), ...)`
- `predict.qregression` -> `eval(npqreg(...), ...)`
- `predict.singleindex` -> `eval(npindex(...), ...)`
- `predict.smoothcoefficient` -> `eval(npscoef(...), ...)`

3. Bandwidth-object `predict.<bwclass>()` methods that invoke estimators:
- `predict.bandwidth` -> `npudens(...)`
- `predict.dbandwidth` -> `npudist(...)`
- `predict.rbandwidth` -> `npreg(...)`
- `predict.conbandwidth` -> `npcdens(...)`
- `predict.condbandwidth` -> `npcdist(...)`
- `predict.plbandwidth` -> `npplreg(...)`
- `predict.sibandwidth` -> `npindex(...)`
- `predict.scbandwidth` -> `npscoef(...)`

4. `se()` / `gradients()` generic paths:
- Generic wrappers are in `R/util.R`; most methods are accessors, but their parent objects are created through the distributed estimator pathways above.

## 4) MPI Runtime Functions Relevant to Guardrails
Strong candidates to use inside helper/preflight:
- `mpi.comm.size`, `mpi.comm.rank`
- `mpi.barrier`
- `mpi.remote.exec` (for slave-side existence/options checks)
- `mpi.bcast.Robj2slave` (auto-export)
- `mpi.bcast.cmd(..., caller.execute=TRUE)` (internal execution primitive)
- `mpi.comm.set.errhandler` (errors-return behavior)
- `mpi.close.Rslaves(force=TRUE)` (recovery)

## 5) Clobber / Misuse Risks to Explicitly Guard
1. Nested broadcast risk:
- user wraps `np*` call in `mpi.bcast.cmd(..., caller.execute=TRUE)` while auto-dispatch is enabled.
- guard needed: broadcast-context sentinel and no recursive broadcast.

2. Parent-frame evaluation risk:
- many `predict.*` methods call `eval(np*(...), envir=parent.frame())`.
- helper must resolve dependencies from caller frame/formula symbols before distributed execution.

3. Slave global assignment collisions:
- `mpi.bcast.Robj2slave` assigns names directly in slave `.GlobalEnv`.
- helper should reserve a private prefix and avoid temporary-name leakage.

4. Options drift:
- `options(np.messages=..., np.tree=..., np.largeh.rel.tol=..., np.disc.upper.rel.tol=...)` can differ by rank if not synchronized.

5. Partial initialization:
- neither `npRmpi.start()` nor `np.mpi.initialize()` has run in a usable context.
- helper should auto-bootstrap where feasible and otherwise fail fast with exact instruction.

## 6) Implementation Checklist (What Must “Play Nicely”)
Use this checklist when wiring the helper:

1. Intercept all top-level `np*` constructors listed in Section 1.
2. Intercept `npplot` and all `plot.*` methods that route to `npplot`.
3. Intercept all `predict.*` methods that route through `eval(np*(...))`.
4. Intercept bandwidth `predict.*` methods that call estimators.
5. Add nested-broadcast sentinel handling.
6. Add preflight + dependency discovery + auto-broadcast + option sync.
7. Ensure graceful recovery path and actionable errors (no silent hang).

## 7) Recommendation on Scope
For first implementation pass, prioritize:
1. `npregbw`, `npreg`, `npplot`, `predict.npregression`, `predict.rbandwidth`
2. `npcdensbw`, `npcdens`, `predict.condensity`, `predict.conbandwidth`
3. `npcdistbw`, `npcdist`, `predict.condistribution`, `predict.condbandwidth`
4. then remaining estimator families.

This staged rollout covers the highest-frequency paths and most observed hang scenarios first.

## 8) Exhaustiveness Audit Addendum (New Findings)
Additional scrutiny uncovered several omission classes that should be in the master list.

### 8.1 Native C Symbol Surface
`/Users/jracine/Development/np-npRmpi/src/np_init.c` registers `.C` entry points and leaves dynamic lookup enabled:

1. Registered `.C` APIs (high-impact for estimator execution):
- `np_density`, `np_density_bw`
- `np_density_conditional`, `np_density_conditional_bw`
- `np_distribution_bw`, `np_distribution_conditional_bw`
- `np_regression`, `np_regression_bw`
- `np_quantile_conditional`
- `np_kernelsum`
- `np_mpi_init`, `np_release_static_buffers`, `np_set_seed`, `np_set_tgauss2`
- `gsl_bspline`, `gsl_bspline_deriv`

2. Dynamic `.Call` lookup is explicitly enabled:
- `R_useDynamicSymbols(dll, TRUE)`
- This keeps legacy `.Call("mpi_*", ...)` symbols callable and broadens the runtime surface.

Implication:
- Guardrail design must include native-call failure handling and initialization-state checks before any path that eventually reaches these symbols.

### 8.2 MPI Runtime Names Not in Initial “np*” List
The exported MPI API is large (100+ names) and can interfere with helper assumptions if users call them directly mid-workflow:
- transport/state: `mpi.send`, `mpi.recv`, `mpi.probe`, `mpi.iprobe`, `mpi.wait*`, `mpi.test*`
- topology/communicator: `mpi.comm.*`, `mpi.intercomm.merge`, `mpi.comm.set.errhandler`
- broadcast/object transfer: `mpi.bcast.*`, `mpi.remote.exec`, `mpi.scatter.*`, `mpi.gather.*`

Implication:
- helper should include a session health preflight before each distributed `np*` call, because user-side low-level MPI calls may desynchronize ranks.

### 8.3 Clobber-Prone Global Temp Symbols (Concrete)
Current R MPI utilities create/expect global temporary names on slaves:
- `.tmpname`, `.tmpRobj`, `.tmp.obj`, `.tname`, `.tinfo`
- `.mpi.err`, `.mpi.applyLB` (historical/active references)

Sources:
- `/Users/jracine/Development/np-npRmpi/R/Rcoll.R`
- `/Users/jracine/Development/np-npRmpi/R/Rparutilities.R`

Implication:
- new helper should reserve a dedicated namespace prefix (e.g. `.__npRmpi_*`) and clean up after use.
- add collision checks before assignment in slave `.GlobalEnv`.

### 8.4 Exported Dot-Prefixed Internal Worker Functions
Several internals are exported and thus user-callable:
- `.mpi.worker.exec`, `.mpi.worker.apply`, `.mpi.worker.sim`, `.mpi.worker.applyLB`
- `.docall`, `.splitIndices`, `.simplify`, `.typeindex`, `.force.type`, `.mpi.undefined`

Implication:
- accidental user calls can perturb worker state.
- ideally keep compatibility but treat these as reserved and document “do not call directly”.

### 8.5 Additional Generic/Method Surface Beyond Initial Focus
The initial list emphasized constructors/plot/predict/se/gradients. Additional relevant pathways include:
- `fitted.*`, `residuals.*`, `coef.*`, `vcov.*`, `summary.*`, `print.*` for many classes
- Most are accessors/formatters, but some trigger estimator re-evaluation transitively via `predict.*`/`plot.*` paths.

Implication:
- dispatcher integration tests should include at least one representative call for each generic family, not only constructors.

### 8.6 Session/Option Side Effects to Track
Package load/startup mutates options:
- `np.messages`, `np.tree`, `np.largeh.rel.tol`, `np.disc.upper.rel.tol`, `np.groupcv.fast`
- `npRmpi.reuse.slaves` defaulting logic

Source:
- `/Users/jracine/Development/np-npRmpi/R/zzz.R`

Implication:
- guardrail helper should sync a controlled whitelist of options to slaves and verify parity.

## 9) Expanded “Must-Not-Clobber” Master Checklist
1. User objects referenced by formula and pass-through arguments.
2. Reserved helper temp symbols (`.__npRmpi_*` preferred; avoid existing `.tmp*` collisions).
3. MPI worker state (`.mpi.err`, pending daemon loop semantics).
4. Global options affecting kernel/CV paths.
5. Existing slave pool lifecycle (`npRmpi.reuse.slaves`, `npRmpi.start/stop` semantics).
6. S3 transitive paths that re-enter constructors (`predict.*`, `plot.*`, bandwidth `predict.*`).

## 10) Recommended Additional Validation Before Implementation
1. Static check: enumerate every `predict.*` and `plot.*` method body and classify as:
- accessor-only
- re-evaluates estimator (`eval(np*...)`)
- calls `npplot`.
2. Runtime check harness:
- for each top-level `np*` constructor + one method chain, run:
  - no manual broadcast (expected success),
  - manual nested `mpi.bcast.cmd` wrap (expected trap/warn),
  - missing object on slave (expected fast error message),
  - options mismatch (expected auto-sync).

## 11) End-to-End Chain Trace (Master -> Slave -> Native)
This section maps key production chains from user entry to native/MPI layers.

### 11.1 Regression bandwidth/fit chain
1. User call:
- `npregbw(...)` in `/Users/jracine/Development/np-npRmpi/R/np.regression.bw.R`
2. Native call:
- `.C("np_regression_bw", ..., PACKAGE="npRmpi")` in `/Users/jracine/Development/np-npRmpi/R/np.regression.bw.R`
3. C entry:
- `void np_regression_bw(...)` in `/Users/jracine/Development/np-npRmpi/src/np.c`
4. Objective/CV internals:
- function pointer assignment by `bwmethod` then calls into CV routines in `/Users/jracine/Development/np-npRmpi/src/jksum.c`
5. MPI collectives:
- reduction/gather calls in `/Users/jracine/Development/np-npRmpi/src/jksum.c` (`MPI_Allreduce`, `MPI_Gather`, `MPI_Bcast`, etc.).
6. Fit call:
- `npreg(...)` in `/Users/jracine/Development/np-npRmpi/R/np.regression.R` -> `.C("np_regression", ...)`
7. Transitive re-entry:
- `predict.npregression()` in `/Users/jracine/Development/np-npRmpi/R/regression.R` calls `eval(npreg(...), envir=parent.frame())`
- `plot.npregression()` calls `npplot(...)`

### 11.2 Conditional density bandwidth/fit chain
1. User call:
- `npcdensbw(...)` in `/Users/jracine/Development/np-npRmpi/R/np.condensity.bw.R`
2. Native call:
- `.C("np_density_conditional_bw", ..., PACKAGE="npRmpi")`
3. C entry:
- `void np_density_conditional_bw(...)` in `/Users/jracine/Development/np-npRmpi/src/np.c`
4. CV internals:
- `cv.ml`/`cv.ls` objective dispatch in `np.c` to `jksum.c` kernels.
5. Fit call:
- `npcdens(...)` in `/Users/jracine/Development/np-npRmpi/R/np.condensity.R` -> `.C("np_density_conditional", ...)`
6. Transitive re-entry:
- `predict.condensity()` in `/Users/jracine/Development/np-npRmpi/R/condensity.R` calls `eval(npcdens(...))`
- `plot.condensity()` calls `npplot(...)`

### 11.3 Conditional distribution bandwidth/fit chain
1. User call:
- `npcdistbw(...)` in `/Users/jracine/Development/np-npRmpi/R/np.condistribution.bw.R`
2. Native call:
- `.C("np_distribution_conditional_bw", ..., PACKAGE="npRmpi")`
3. C entry:
- `void np_distribution_conditional_bw(...)` in `/Users/jracine/Development/np-npRmpi/src/np.c`
4. CV internals:
- objective dispatch to distribution CV routines in `jksum.c`.
5. Fit call:
- `npcdist(...)` in `/Users/jracine/Development/np-npRmpi/R/np.condistribution.R` uses `.C("np_density_conditional", ..., densOrDist=NP_DO_DIST, ...)`
6. Transitive re-entry:
- `predict.condistribution()` in `/Users/jracine/Development/np-npRmpi/R/condistribution.R` calls `eval(npcdist(...))`
- `plot.condistribution()` calls `npplot(...)`

### 11.4 Slave execution loop chain
1. Startup/init:
- `npRmpi.start()` in `/Users/jracine/Development/np-npRmpi/R/session.R` calls `mpi.spawn.Rslaves(...)` and `mpi.bcast.cmd(np.mpi.initialize(), caller.execute=TRUE)`
2. Slave daemon:
- `/Users/jracine/Development/np-npRmpi/inst/slavedaemon.R` loops:
  - `tmp.message <- mpi.bcast.cmd(...)`
  - `eval(tmp.message, envir=.GlobalEnv)`
3. Broadcast primitive:
- `mpi.bcast.cmd()` in `/Users/jracine/Development/np-npRmpi/R/Rcoll.R` serializes command+args, sends to each rank, optionally executes on caller.

## 12) Newly Identified Gaps / Misses
1. Daemon-context recursion trap still required:
- If auto-dispatch wraps code already executing inside `mpi.bcast.cmd` on slaves, nested broadcast can deadlock; explicit sentinel detection is required in both master and slave contexts.

2. Health-check must include communicator/rank alignment, not just "MPI initialized":
- `np_mpi_init` records rank/size, but user low-level MPI calls can desynchronize communicator state after initialization.

3. `predict.*` and bandwidth `predict.*` are mandatory interception points:
- They are common hidden re-entry points (`eval(np*...)`) and will bypass a constructor-only dispatcher.

4. `plot.* -> npplot -> np*` needs one-shot dispatch semantics:
- `npplot` can loop across slices and call `np*` repeatedly; helper must avoid rebroadcasting each nested step.

5. Option synchronization must be explicit and whitelisted:
- `np.messages`, `np.tree`, `np.largeh.rel.tol`, `np.disc.upper.rel.tol`, `np.groupcv.fast` affect kernel/CV branches and must be synchronized before estimator entry.

6. Temporary symbol namespace policy is still missing:
- Existing `.tmp*` globals in MPI utility code increase collision risk; helper temp objects should use a dedicated prefix and cleanup lifecycle.

## 13) Cross-Repo Audit (np-master vs np-npRmpi)
Cross-reference was run to catch callable surfaces potentially missed by focusing only on `np-npRmpi`.

### 13.1 Exported top-level API parity
Result:
1. `np-master` exported API is a strict subset of `np-npRmpi`.
2. No non-MPI top-level `np*` exports were missing from `np-npRmpi`.
3. Extra exports in `np-npRmpi` are MPI/runtime utilities only.

Implication:
- our constructor-level inventory is complete for top-level `np*` functions.

### 13.2 S3 registration differences that matter
`np-master` has these S3 registrations that `np-npRmpi` does not register in `NAMESPACE`:
1. `print.npregiv`
2. `summary.npregiv`
3. `plot.npregiv`
4. `print.npregivderiv`
5. `summary.npregivderiv`
6. `plot.npregivderiv`
7. `gsl.bs.default`
8. `predict.gsl.bs`

Notes:
1. The `npregiv*` method functions exist in `np-npRmpi/R/npregiv.R` and `np-npRmpi/R/npregivderiv.R`, but are not S3-registered in `np-npRmpi/NAMESPACE`.
2. `gsl.bs` / `predict.gsl.bs` functions exist in both repos, but only `np-master` registers the related S3 methods.

Implication for helper/guardrail scope:
1. Add `npregiv` and `npregivderiv` constructor chains to explicit coverage tests (already exported and example-used).
2. Even if their `plot.*` methods are base-plot only (not `npplot`), include them in callable inventory so user-exposed behavior is not omitted.
3. Track the S3-registration drift as a separate namespace-consistency issue.

### 13.3 Example-surface cross-check
`man/*.Rd` example parsing confirms both repos exercise the same `np*` family, including:
1. `npregiv`
2. `npregivderiv`
3. all major bandwidth/fit/test functions listed in Section 1.

Implication:
- inventory should explicitly mark `npregiv`/`npregivderiv` as example-critical callable chains.

## 14) Current Coverage Status (Auto-Dispatch)
Implemented constructor/method interception currently includes:
1. `npudensbw`, `npudens`
2. `npudistbw`, `npudist`
3. `npregbw`, `npreg`
4. `npcdensbw`, `npcdens`
5. `npcdistbw`, `npcdist`
6. `npscoefbw`, `npscoef`
7. `npindexbw`, `npindex`
8. `npplregbw`, `npplreg`
9. `npconmode`
10. `npksum`
11. `npcopula`
12. `npquantile`
13. `np.pairs`, `np.pairs.plot`
14. `npplot`
15. `npregiv`, `npregivderiv`
16. `npqreg`
17. `npcmstest`, `npqcmstest`, `npdeneqtest`
18. `npdeptest`, `npsdeptest`, `npsymtest`, `npunitest`, `npsigtest`

Still pending for equivalent treatment:
1. confirm end-to-end method chains for accessor-driven re-entry (`predict.*`/`plot.*`) beyond currently wrapped constructors.
2. validate desired semantics for unconditional helper utilities (`npseed`, `nptgauss`, `uocquantile`, `npuniden.*`) that are not distributed estimators but can appear in user pipelines.

## 15) Known Open Issue
1. `npplreg` crash mode observed when called without an active slave pool (`npRmpi.start(...)` not invoked), with native trace reaching `kernel_weighted_sum_np` in `jksum.c`.
2. Early guardrail added in R-layer for `npplregbw*` and `npplreg*` to fail fast with a startup instruction instead of entering undefined native execution.
3. Keep this tracked separately as `npplreg crash issue` for focused C-level forensic analysis; do not block broader user-facing dispatch/guardrail rollout.

## 16) Transitive Audit Status (predict/plot chains)
Static audit status:
1. All `predict.*` and bandwidth `predict.*` methods that re-enter estimators (`eval(np*...)`) now route through constructor entry points with both:
- `.npRmpi_require_active_slave_pool(...)`
- `.npRmpi_autodispatch_call(...)` (when active)
2. All `plot.*` methods that delegate to `npplot(...)` are covered by `npplot()` guard + autodispatch.
3. No remaining autodispatch call sites were found without a nearby active-pool guard (excluding helper internals in `R/np.autodispatch.R` by design).
