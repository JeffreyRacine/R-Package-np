# Performance Benchmarks

Use this folder for timing-focused scripts and summary utilities.

- `methods/`: estimator-family benchmark suites (`npreg`, `npcdens`, etc.).
- `oneoff/`: one-off benchmark harnesses for methods not covered by combo grids.

## Output Conventions

- Write raw and summary artifacts to `/tmp`.
- Include key timing/objective columns (`elapsed`, objective value where relevant,
  selected bandwidths, and `num.feval` when available).
- For post-vs-pre reports, use positive percent values to mean "post is faster".
