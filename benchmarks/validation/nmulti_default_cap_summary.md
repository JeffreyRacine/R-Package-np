# `nmulti` Default-Cap Decision

These validation harnesses support the default multistart-cap change
from `min(5, p)` to `min(2, p)`.

Study components in this folder:

- `nmulti_core_surface_study.R`
  - synthetic `npreg`/`npcdens` surface-stability study spanning
    `lc`, `ll`, and `lp` with `degree=1,2`
- `nmulti_borderline_density_dist_study.R`
  - targeted `npcdens`/`npcdist` follow-up on the previously borderline
    `ll` and `lp(degree=2)` routes under mixed-data DGPs
- `nmulti_wage1_start_stress.R`
  - real-data `wage1` LP regression stress test over varied starting
    bandwidth vectors
- `nmulti_wage1_nomad_stress.R`
  - real-data `wage1` LP automatic-degree stress test for
    `search.engine="nomad+powell"`

High-level conclusion:

- `nmulti=5` did not show enough empirical benefit to justify keeping
  the older default cap.
- `nmulti=1` was often adequate on synthetic designs, but it showed
  clear bad-tail brittleness on real-data LP and NOMAD start-stress
  problems.
- `nmulti=2` consistently tracked the practical `nmulti=5` reference
  far more closely while still yielding large runtime savings.

Recommended default:

- use `min(2, p)` for the core bandwidth selectors
- use `nmulti=1` for the univariate boundary density helper

Raw reports and CSV summaries are written to dated `/tmp` artifact roots
when the scripts are run.
