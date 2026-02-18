# AGENTS.md

Use `/Users/jracine/Development/AGENTS.md` as the canonical workflow and benchmarking policy.

Repo-specific note:
- Merge (do not overwrite) R-layer and man-page files when porting from `np-master`.
- Any new or ported benchmark/function must follow MPI protocol from `/Users/jracine/Development/AGENTS.md`: start/stop slaves, broadcast commands/options/data explicitly, support passable `nslaves`/`rslaves`, and keep parity with serial seed/DGP/options.
- Keep shared regression summary telemetry in sync with `np-master`: `num.feval`, `num.feval.fast`, `num.feval.fallback`.
- Keep `man/np.kernels.Rd` and kernel/options/plot cross-links aligned, while preserving package-specific `npRmpi` docs text.
