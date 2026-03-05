# npRmpi Native Registration Policy Decision (2026-03-01)

## Scope
- Package: `npRmpi`
- Files reviewed:
  - `NAMESPACE`
  - `src/np_init.c`
  - `R/Rmpi.R`, `R/Rcomm.R`, `R/Rcoll.R` (legacy MPI `.Call("mpi_*", ...)` entrypoints)

## Current Facts
1. `NAMESPACE` uses `useDynLib(npRmpi)` (no `.registration=TRUE` flag).
2. `src/np_init.c` registers `C_*` routines with `R_registerRoutines(...)`.
3. `src/np_init.c` currently calls `R_useDynamicSymbols(dll, TRUE)`.
4. The R layer still uses many string-based `.Call("mpi_*", PACKAGE="npRmpi")` entrypoints inherited from embedded Rmpi APIs.

## Decision
Keep the current native symbol policy for this release cycle:
- keep `useDynLib(npRmpi)` as-is,
- keep `R_useDynamicSymbols(dll, TRUE)` as-is,
- do not flip to strict native-registration mode in this tranche.

## Rationale
1. Flipping to strict registration now would break unresolved legacy MPI `.Call("mpi_*", ...)` calls unless all such symbols are explicitly registered and validated first.
2. That migration is a larger compatibility change and not a housekeeping-only edit.
3. This tranche is constrained to no-regression housekeeping/validation work.

## Follow-up Plan (separate tranche)
1. Build complete inventory of all `mpi_*` native symbols used by R wrappers.
2. Add explicit registration entries for required `mpi_*` routines.
3. Update `NAMESPACE` to `useDynLib(npRmpi, .registration=TRUE)`.
4. Flip `R_useDynamicSymbols(dll, FALSE)` only after full route validation (session/attach/manual/profile) and `R CMD check --as-cran` pass.

## Validation Notes
- No runtime behavior change introduced by this policy decision document.
- Existing check baseline remains unchanged (`checkbashisms` WARNING + CRAN incoming version NOTE).

## Canonical Implementation Directive (2026-03-05)

This repository follows a strict canonical execution rule:

1. One canonical implementation per method (outside explicit `np.tree` branching).
2. Unsupported configurations must fail fast with explicit `stop(...)` diagnostics.
3. No silent remap/coercion of user-selected options (for example `bwmethod`, `regtype`, kernels, `cv.iterate`, or bounds transforms).
4. No hidden alternate execution paths for the same method semantics.
5. All fit-defining options (for example `degree`, `basis`, `bernstein.basis`, kernels, and bounds) must be propagated and used by the canonical path.
6. `np.tree=FALSE` is the default; when `np.tree=TRUE`, behavior must remain semantics-preserving and option-compatible with the canonical path.
7. Remove or reject legacy/debug compatibility branches that add redundant runtime overhead once canonical behavior exists.
