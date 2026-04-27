# Scale-Factor Floor Integrated Contract Plan

Date: 2026-04-26

Scope:

- Primary repo: `/Users/jracine/Development/np-master`
- `npRmpi`: read-only until `np-master` has a diagnosed, installed-build-gated contract
- No manuscript or simulation campaign files

Approved decision record:

- `/Users/jracine/Development/tmp/scale_factor_floor_integrated_contract_20260426/DECISION_RECORDS_20260426.md`

## Goal

Create a complete, consistent, scientifically grounded scale-factor and floor
contract across bandwidth-selection routes, with the highest chance of fixing
the real problem and the lowest chance of regression or collateral damage.

The target is not merely to patch the latest `npcdens` symptom. The target is to
make initialization controls, optimizer feasibility, final acceptance, returned
objects, and documentation all agree with estimator-specific bandwidth scaling.

## Critique Of The Prior Plan

The previous plan was directionally right but not yet strong enough for this
risk level.

1. It moved too quickly toward shared helpers.
   - A shared helper can reduce duplication, but only after we prove the route
     contracts are actually shared. Premature consolidation could spread one
     mistaken interpretation across every estimator family.
2. It did not make the audit itself a hard gate.
   - The plan should prohibit code changes beyond the current `npcdens`
     candidate until a route-by-route ledger identifies which paths are
     structurally equivalent and which are not.
3. It under-specified the scaling source of truth.
   - It stated that hard-coded exponents are forbidden, but did not require a
     concrete inventory of where each family currently computes `scale_j` and
     `nconfac`, or how those values cross the R/C boundary.
4. It blurred "floor-aware starts" with "floor-constrained search".
   - Making `cfac.init` floor-aware may be correct, but it must be justified as
     part of an initialization contract, not as a substitute for optimizer
     bounds or final feasibility.
5. It did not separate eval-only, explicit-bandwidth, normal-reference, Powell,
   NOMAD, and NOMAD+Powell enough.
   - These routes have different obligations. A floor used for search may not
     be appropriate for storing or evaluating explicit user bandwidths.
6. It did not define "invalid objective" precisely enough.
   - We need a raw-objective versus wrapped-objective ledger. Penalties may guide
     the optimizer, but only raw objective validity can certify a final optimum.
7. It did not include enough negative tests.
   - We need tests for bad `hbc.init`, below-floor starts, infeasible explicit
     floors, and transformed/untransformed scale mismatches, not only successful
     high-floor cases.
8. It made documentation a late phase.
   - Documentation wording should be drafted before broad implementation so the
     code is held to a human-readable contract.
9. It did not specify rollback criteria per tranche.
   - Each candidate patch must have a predeclared rollback trigger.
10. It did not explicitly preserve current user-facing behavior where the floor
    is not active.
    - Baseline-equivalent runs with `scale.factor.lower.bound = 0` or a
      permissive floor must be part of every family gate.

## Scientific Contract

For each continuous bandwidth component:

```text
h_j = sf_j * scale_j * nconfac
```

where:

- `sf_j` is the dimensionless scale factor being initialized, searched, and
  reported when the route is in scale-factor mode;
- `scale_j` is the route-appropriate data scale, using the package's established
  convention, typically `EssDee` for the relevant continuous variable;
- `nconfac` is estimator-specific and must be computed from the actual kernel
  order and the relevant continuous dimension for that estimator family;
- tests and helpers must derive `nconfac` from object metadata or existing route
  helpers, not from hard-coded exponents such as `n^(-1/6)`.

The contract separates five concepts:

1. `cfac.init`
   - deterministic first continuous scale factor used for the first search
     start.
2. `lbc.init`
   - lower endpoint for random continuous scale-factor draws used in multistart
     initialization.
3. `hbc.init`
   - upper endpoint for random continuous scale-factor draws used in multistart
     initialization.
4. `scale.factor.lower.bound`
   - minimum admissible continuous scale factor during fixed-bandwidth search.
5. raw objective validity
   - the final selected candidate must have a finite raw objective, independent
     of any wrapped penalty used for optimization guidance.

## Approved Decisions

Jeffrey approved the following contract decisions in the 2026-04-26 thread.
They govern the implementation plan unless superseded by a later explicit
decision.

1. Inconsistent random-start ranges error clearly.
   - If `hbc.init < max(lbc.init, scale.factor.lower.bound)`, the route must
     fail clearly rather than silently widening `hbc.init`.
2. Deterministic first search starts are floor-aware.
   - The effective deterministic first continuous scale factor is
     `max(cfac.init, scale.factor.lower.bound)`.
3. Explicit user bandwidths remain user-controlled outside search starts.
   - A user may provide any explicit bandwidth when it is not being used as a
     starting value for search.
4. Invalid raw objectives fail.
   - A selected final candidate with an invalid raw objective must produce a
     clear error, not a successful object with a penalty/sentinel objective.

## High-Standard Operating Rules

1. Evidence precedes abstraction.
   - Do not introduce or broaden shared helpers until the audit ledger proves
     that affected routes share the same contract.
2. One risk axis per tranche.
   - Do not combine validation, optimizer feasibility, final acceptance,
     documentation, and cross-family refactor in one patch.
3. No silent semantics changes.
   - If a user-supplied control is inconsistent with the floor, fail clearly
     unless Jeffrey explicitly chooses a widening/remediation policy.
4. No final penalty optima.
   - Wrapped penalties can guide Powell/NOMAD but cannot certify a successful
     bandwidth object.
5. No hard-coded scaling shortcuts.
   - Scaling must come from route metadata, existing helper state, or a new
     audited metadata helper.
6. Preserve eval-only semantics unless explicitly changed.
   - Explicit user bandwidth storage/evaluation is not automatically the same as
     search.
7. Installed-build proof governs claims.
   - Source-loaded probes are diagnostic only.
8. `np-master` first.
   - `npRmpi` remains read-only until `np-master` passes gates.

## Additional Regression-Prevention Controls

These controls are deliberately conservative. They are intended to keep the
work boring, auditable, and reversible.

1. Use a change ladder.
   - Level 0: read-only audit and scripts.
   - Level 1: tests/sentinels only.
   - Level 2: one-route code change.
   - Level 3: shared helper introduced but used by one route only.
   - Level 4: one additional route adopts the shared helper.
   - Level 5: broader adoption after two independent route families pass.
2. Require a stop-the-line trigger.
   - Any contradiction from a simple repro stops implementation and returns the
     work to diagnosis.
3. Keep a pre/post golden-master layer.
   - For every touched route, compare permissive-floor behavior against the
     pre-floor or baseline-equivalent installed build wherever admissible.
4. Use decision records.
   - Every semantic choice gets a short dated decision record under the audit
     root before code changes.
5. Prefer additive checks before rewrites.
   - Add diagnostics, tests, and narrow guards before replacing architecture.
6. Keep package-code patches smaller than the audit.
   - If a patch is large enough that it cannot be explained by one ledger row
     and one failure mechanism, split it.
7. Make error behavior part of the public contract.
   - A clear intentional error is acceptable only when it is documented,
     tested, and less misleading than returning an invalid object.
8. Keep simulation/campaign evidence out of the acceptance path.
   - Campaign sentinels are useful regression probes, but package acceptance
     must rest on package-local installed-build gates.

## Phase 0: Freeze And Protect

Before more code changes:

1. Preserve the current `npcdens` forensic bundle:
   - `/Users/jracine/Development/tmp/npcdens_scale_floor_repair_20260426`
2. Record current working-tree edits and mark them as a candidate, not a final
   integrated solution.
3. Create a new audit root:
   - `/Users/jracine/Development/tmp/scale_factor_floor_integrated_contract_20260426`
4. Do not touch `.Rd` files, other families, `npRmpi`, manuscript files, or
   simulation campaign files in this phase.
5. Save a patch snapshot of the current `npcdens` candidate so it can be
   reverted or split without losing work.
6. Produce an explicit "current dirty tree is not rebuild-safe" note if
   uncommitted candidate code remains in `np-master`.

Exit gate:

- a short status note states exactly what is already candidate-patched and what
  remains unaudited.
- the current candidate patch is saved under the audit root.

## Phase 1: Build The Contract Ledger

Create a route-by-route ledger before patching.

Rows:

- `npregbw`
- `npudensbw`
- `npcdensbw`
- `npudistbw`
- `npcdistbw`
- `npscoefbw`
- `npindexbw`
- `npplregbw`
- internal hot-start routes
- eval-only routes
- normal-reference routes
- NOMAD-only routes
- NOMAD+Powell handoff routes

Columns:

1. public initialization controls exposed;
2. default values;
3. current R validation for `cfac.init`;
4. current R validation for `lbc.init`;
5. current R validation for `hbc.init`;
6. current R validation for `scale.factor.lower.bound`;
7. current policy when `hbc.init < lbc.init`;
8. current policy when `hbc.init < max(lbc.init, floor)`;
9. R-side effective deterministic start;
10. R-side effective random lower;
11. R-to-C payload position and name;
12. matching C enum/index;
13. C initialization function used;
14. scale convention for each continuous variable;
15. `nconfac` formula source;
16. transformed/untransformed search state;
17. objective wrapper validity checks;
18. final acceptance checks;
19. returned object scale;
20. explicit eval-only behavior;
21. current `.Rd` wording;
22. current sentinel coverage.

Ledger requirements:

- every formula cell must cite the exact file and line where the value is
  computed or passed;
- unknown cells must be marked `UNKNOWN`, not inferred;
- structurally distinct routes must not be grouped prematurely.
- build the ledger from a small script where possible, then manually review it;
- the ledger must include both R-layer and C-layer evidence;
- the ledger must classify each route as:
  - `shared contract likely`,
  - `shared contract possible with exception`,
  - `route-specific contract`,
  - `unknown`.

Exit gate:

- no package code changes are allowed until the ledger is complete for
  `npcdens`, `npregbw`, `npudensbw`, `npudistbw`, `npcdistbw`, and at least
  read-only mapped for `npindexbw`, `npscoefbw`, and `npplregbw`.
- no route may be patched from `UNKNOWN`.

## Phase 1.5: Build Golden-Master Baselines

Before changing any additional route, build installed golden-master outputs.

For each route with a sentinel:

- baseline-equivalent or permissive floor;
- default floor;
- high floor;
- intentionally bad `hbc.init`;
- eval-only explicit below-floor bandwidth where applicable.

For each result record:

- package path;
- package version;
- source commit or tarball SHA;
- `.libPaths()`;
- call;
- seed/data identifier;
- raw objective if available;
- wrapped objective if available;
- returned objective;
- bandwidth vector;
- scale-factor vector where available;
- floor vector derived from route metadata;
- validity status.

Exit gate:

- the baseline table is saved before any non-`npcdens` patch.
- any expected change in behavior is predeclared before the candidate is run.

## Phase 2: Define The Public Contract In Words

Before implementation, draft user-facing wording for all affected `*bw.Rd`
files.

Required definitions:

- `cfac.init`: deterministic first continuous scale factor used to form the
  first search bandwidth. It is an initialization control, not the admissibility
  floor.
- `lbc.init`: lower endpoint for random continuous scale-factor draws used by
  multistart initialization.
- `hbc.init`: upper endpoint for random continuous scale-factor draws used by
  multistart initialization.
- `scale.factor.lower.bound`: admissibility floor for searched continuous scale
  factors in fixed-bandwidth search.

Required clarifications:

- random starts use an effective lower endpoint no smaller than the active
  floor;
- deterministic starts may be adjusted to be feasible under the active floor,
  but that does not make `cfac.init` the floor;
- all physical bandwidth floors use estimator-specific scaling and continuous
  dimension;
- explicit eval-only bandwidth behavior is documented separately from search.

Exit gate:

- wording is reviewed against the ledger before any `.Rd` edits are made.

## Phase 3: Implement Approved Inconsistent-Range Policy

Implement this before shared validation helpers:

Case:

```text
hbc.init < max(lbc.init, scale.factor.lower.bound)
```

Approved policy:

- fail clearly with an error naming `hbc.init`, `lbc.init`,
  `scale.factor.lower.bound`, and the effective lower endpoint.

Rejected default policy:

- silently widening `hbc.init`, because it changes the user-selected random
  multistart distribution.

Possible exception:

- internal hot-start routes may use non-public zero-valued placeholders if the
  route is not generating random starts. Such exceptions must be named in the
  ledger and tested.

Exit gate:

- the written decision record is retained under the audit root;
- tests cover the approved error policy and any named exceptions.

## Phase 4: Add Diagnostic Sentinels Before More Repair

Add or retain scratch sentinels before broad package changes.

Sentinel classes:

1. baseline-equivalent permissive-floor sentinel;
2. default-floor sentinel;
3. high-floor sentinel;
4. invalid raw-objective sentinel;
5. bad `hbc.init` sentinel;
6. transformed-bound sentinel;
7. eval-only explicit below-floor bandwidth sentinel;
8. NOMAD-only sentinel;
9. NOMAD+Powell handoff sentinel.

Rules:

- each family gets one small fast sentinel before any family patch;
- `npcdens` keeps the heavier patent and production sentinels;
- raw logs and CSVs are retained under the audit root;
- sentinel scripts must record package path and `.libPaths()`.
- sentinels must be executable against both baseline and candidate installed
  libraries without editing the script;
- tests must include at least one negative assertion per changed contract.

Exit gate:

- every planned patch has a named pre/post sentinel.

## Phase 5: Stabilize `npcdens` Candidate

Use `npcdens` as the reference implementation only after its remaining edge
contracts are verified.

Required checks:

- no hard-coded exponent remains in tests or diagnostics where object metadata
  is available;
- `cfac.init`, `lbc.init`, `hbc.init`, and `scale.factor.lower.bound` are
  validated coherently;
- bad `hbc.init` range fails clearly;
- Powell rejects below-floor candidates during search;
- final acceptance recomputes and verifies raw objective;
- `scale.factor.lower.bound = 0` or a permissive floor recovers admissible
  baseline behavior;
- eval-only explicit bandwidth behavior is unchanged and documented.

Rollback triggers:

- patent permissive-floor LC no longer matches the Apr 9 admissible optimum;
- bounded production sentinels regress;
- eval-only explicit bandwidth tests change behavior unexpectedly;
- error text becomes vague or masks the failing condition.

Exit gate:

- installed `npcdens` sentinel set passes.
- a written decision states whether the current `npcdens` candidate remains one
  tranche or must be split into:
  - initialization validation,
  - objective-wrapper feasibility,
  - final raw-objective acceptance.

## Phase 6: Patch Only Proven Shared Validation

Only after the ledger proves commonality, introduce small R helpers.

Candidate helper responsibilities:

- validate scalar finite controls;
- compute effective deterministic start;
- compute effective random lower;
- validate `hbc.init` against the effective lower;
- return a named list for explicit R-to-C payload construction.

Non-responsibilities:

- do not compute physical floors unless the helper has route metadata;
- do not decide eval-only behavior;
- do not silently modify user-selected `hbc.init`;
- do not alter categorical/count bandwidth rules.

Exit gate:

- helper unit tests pass independently;
- no family is switched to the helper until that family's sentinel set exists.
- helper adoption is limited to one route until that route passes installed
  gates.

## Phase 7: Family-by-Family Adoption

Adopt the contract in this order, one tranche at a time:

1. `npregbw`
2. `npudensbw`
3. `npudistbw`
4. `npcdistbw`
5. `npscoefbw`
6. `npindexbw`
7. `npplregbw`

For each family:

- compare against the ledger;
- decide whether the shared helper applies;
- patch only that family;
- run its sentinel set from an installed build;
- run adjacent existing tests;
- update only that family's documentation after code behavior is gated.
- compare permissive-floor outputs against the golden-master baseline;
- compare default-floor outputs against expected behavior declared before the
  patch;
- record whether any change is intended, incidental but acceptable, or a
  regression.

Rollback triggers:

- any family sentinel changes outside the predeclared expectation;
- explicit bandwidth behavior changes without a documented contract decision;
- a shared helper requires a family-specific exception not listed in the ledger.
- any newly introduced warning/error in an untouched route.

## Phase 8: Documentation Updates

After a family passes installed gates, update its `.Rd` text.

Documentation must:

- distinguish starts from floors;
- define random multistart use of `lbc.init` and `hbc.init`;
- define deterministic use of `cfac.init`;
- define floor enforcement during search and final acceptance;
- specify whether eval-only explicit bandwidths are clamped, rejected, or stored
  unchanged;
- avoid claiming a universal exponent.

Exit gate:

- examples and usage remain unchanged unless explicitly intended;
- `R CMD check` documentation checks pass.

## Phase 9: Final `np-master` Gates

Final gates:

- complete ledger retained under `/Users/jracine/Development/tmp`;
- all family sentinel CSVs retained;
- decision records retained under the audit root;
- candidate patch snapshots retained;
- focused `testthat` files pass from installed build;
- family default/permissive/high-floor sentinels pass;
- `R CMD build`;
- `R CMD check --as-cran`;
- no `ERROR`;
- warnings/notes classified as pre-existing, environmental, or introduced.

Claim discipline:

- call it `candidate keep` after focused installed gates;
- call it `validated` only after the final installed sentinel matrix and
  `R CMD check --as-cran`;
- call it `release-ready` only after tarball-first release checks are accepted.

## Phase 10: `npRmpi` Surgical Parity Port

Only after `np-master` is validated, green-lit, and committed:

- port the accepted `np-master` changes surgically into
  `/Users/jracine/Development/np-npRmpi`;
- do not copy whole files;
- preserve functional equivalence with `np` except for explicit, documented
  MPI-specific execution differences;
- preserve `npRmpi` runtime independence;
- preserve repo-specific R/man behavior;
- update `npRmpi` documentation only where required to keep the public contract
  equivalent;
- run serial and MPI-aware sentinels where applicable;
- verify no runtime `np::` bridge dependency;
- verify no silent serial fallback when MPI execution is selected;
- run `npRmpi` installed-build gates before any parity or release claim.

`npRmpi` exit gate:

- surgical patch diff reviewed against the accepted `np-master` commit;
- `nslaves=1` sentinel behavior matches `np` within the declared numerical
  tolerance for shared routes;
- `nslaves>1` sentinels pass for touched MPI-aware routes;
- final `npRmpi` installed-build checks pass;
- any intentional divergence from `np` is named in a decision record before
  closeout.

## Acceptance Criteria

The integrated contract is acceptable only if:

1. every touched route has a ledger entry with scaling source, `nconfac` source,
   and R/C payload mapping;
2. no helper or test hard-codes estimator exponents when route metadata is
   available;
3. `cfac.init`, `lbc.init`, `hbc.init`, and `scale.factor.lower.bound` are
   validated consistently where the route contract is shared;
4. route-specific exceptions are named and tested;
5. inconsistent random-start ranges fail clearly unless explicitly accepted
   otherwise;
6. below-floor candidates are invalid during search;
7. invalid wrapped penalties cannot become final successful optima;
8. eval-only explicit bandwidth behavior is preserved or deliberately changed
   with documentation and tests;
9. `.Rd` wording matches actual behavior;
10. `np-master` installed-build gates pass before any port or `np` release
    claim;
11. after `np-master` is committed, the accepted changes are surgically ported
    to `npRmpi` and independently validated before closeout;
12. `np`/`npRmpi` functional equivalence is verified for shared routes, with
    any MPI-specific exception explicitly documented.

## Immediate Next Step

Do not broaden the current candidate patch yet.

Next action:

1. produce the contract ledger for `npcdens`, `npregbw`, `npudensbw`,
   `npudistbw`, and `npcdistbw`;
2. identify which routes already have `hbc.init` validation and which do not;
3. apply the approved inconsistent-range decision to the ledger and sentinel
   design before code changes;
4. only then decide whether the current `npcdens` candidate should be kept as
   is, revised, or split into smaller tranches.
