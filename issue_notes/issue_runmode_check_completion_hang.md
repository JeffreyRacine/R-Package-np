# Run-mode `R CMD check --as-cran` completion instability

Observed repeatedly during `man/run` mode checks in `npRmpi`:

- `npRmpi-Ex.Rout` reaches footer (`### * <FOOTER> ... quit('no')`) indicating examples completed.
- `00check.log` can remain truncated at `* checking examples ...` and the check driver may terminate unexpectedly.

This appears related to a known interaction between `Rmpi` and `npRmpi` in check-time process management, not a specific single example failure once current doc patches are applied.

TODO:

1. Reproduce with a minimal package/examples harness outside `npRmpi`.
2. Compare behavior with and without slave reuse (`npRmpi.reuse.slaves`).
3. Compare `R CMD check` with `--no-examples` + standalone examples run to isolate teardown stage.
4. Capture and inspect master/slave process exit paths during check finalization.

## Canonical Implementation Directive (2026-03-05)

This repository follows a strict canonical execution rule:

1. One canonical implementation per method (outside explicit `np.tree` branching).
2. Unsupported configurations must fail fast with explicit `stop(...)` diagnostics.
3. No silent remap/coercion of user-selected options (for example `bwmethod`, `regtype`, kernels, `cv.iterate`, or bounds transforms).
4. No hidden alternate execution paths for the same method semantics.
5. All fit-defining options (for example `degree`, `basis`, `bernstein.basis`, kernels, and bounds) must be propagated and used by the canonical path.
6. `np.tree=FALSE` is the default; when `np.tree=TRUE`, behavior must remain semantics-preserving and option-compatible with the canonical path.
7. Remove or reject legacy/debug compatibility branches that add redundant runtime overhead once canonical behavior exists.
