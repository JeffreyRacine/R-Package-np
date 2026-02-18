# AGENTS.md

Use `/Users/jracine/Development/AGENTS.md` as the canonical workflow and benchmarking policy.

Repo-specific note:
- `np-master` is the source of truth for shared C-core performance changes before porting to `np-npRmpi`.
- If asked for an `npRmpi` counterpart, first verify serial behavior in `np-master` and then port with MPI-specific broadcast/start-stop rules from `/Users/jracine/Development/AGENTS.md` and `/Users/jracine/Development/bench_nprmpi.R`.
- For regression bandwidth summary work, keep `num.feval`, `num.feval.fast`, and `num.feval.fallback` semantics aligned with `np-npRmpi`.
- Keep kernel discoverability docs aligned via `man/np.kernels.Rd` and cross-links in method `.Rd` pages.
