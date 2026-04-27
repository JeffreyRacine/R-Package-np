# Initialization-Parameter Catalogue

Date: 2026-04-27

## Purpose

Catalogue the other `*.init` and initialization-like controls that appear near
the continuous bandwidth-search naming discussion, and decide whether they
should be renamed with `cfac.init`, `lbc.init`, `hbc.init`, and
`scale.factor.lower.bound`.

Conclusion: the continuous scale-factor controls are the only appropriate
renaming target for the current naming tranche. The other initialization-like
controls belong to separate contracts and should remain separate and untouched.

## Search Scope

Repos inspected:

- `/Users/jracine/Development/np-master`
- `/Users/jracine/Development/np-npRmpi`

Search roots:

- `R`
- `man`
- `tests`
- `DESCRIPTION`
- `NAMESPACE`

Search terms:

- public and internal `*.init` names
- `scale.init.categorical.sample`
- `initc.dir`
- `initd.dir`
- `bw.init`
- `npRmpi.init`

## Catalogue And Recommendation

### Continuous Scale-Factor Search Controls

Names:

- `cfac.init`
- `lbc.init`
- `hbc.init`
- `scale.factor.lower.bound`

Role:

- `cfac.init` is the deterministic first continuous scale-factor start.
- `lbc.init` is the lower endpoint for random continuous scale-factor starts.
- `hbc.init` is the upper endpoint for random continuous scale-factor starts.
- `scale.factor.lower.bound` is the admissibility floor for continuous
  fixed-bandwidth search candidates.

Recommendation:

- Eligible for the current naming tranche.
- Rename, if desired, as one coherent continuous scale-factor search family.
- Keep semantics unchanged and route all public names through a shared resolver.

Candidate canonical family:

| Current name | Candidate canonical name |
| --- | --- |
| `cfac.init` | `scale.factor.init` |
| `lbc.init` | `scale.factor.init.lower` |
| `hbc.init` | `scale.factor.init.upper` |
| `scale.factor.lower.bound` | `scale.factor.search.lower` |

### Categorical Scale-Factor Start Controls

Names:

- `dfac.init`
- `lbd.init`
- `hbd.init`

Role:

- categorical analogues for starting categorical bandwidth/lambda search values.
- used in density, distribution, regression, conditional density/distribution,
  and smooth-coefficient routes.
- in smooth-coefficient bandwidth search, these controls have their own
  categorical-specific validity contract, including an upper bound of `2`.

Recommendation:

- Keep separate and untouched in the continuous-search naming tranche.
- Do not force these into the continuous `scale.factor.*` family.
- If renamed later, do it as an independent categorical-search tranche after a
  separate audit of categorical admissible ranges, direction-set interaction,
  and documentation wording.

Possible future names, if a later categorical tranche is explicitly opened:

| Current name | Possible future name |
| --- | --- |
| `dfac.init` | `categorical.factor.init` |
| `lbd.init` | `categorical.factor.init.lower` |
| `hbd.init` | `categorical.factor.init.upper` |

Do not make this change in the continuous-floor naming tranche.

### Powell Direction-Set Initialization Controls

Names:

- `initc.dir`
- `initd.dir`

Role:

- non-random initial values for Powell's direction-set search, separated by
  continuous and categorical variables.
- appear alongside `.dir` controls such as `lbc.dir`, `cfac.dir`, `lbd.dir`,
  `hbd.dir`, and `dfac.dir`.
- exposed in bandwidth selectors and quantile-regression documentation.

Recommendation:

- Keep separate and untouched.
- These are optimizer direction-set controls, not bandwidth or scale-factor
  starting values.
- At most, clarify `.Rd` wording in a documentation-only tranche.

Reason:

- Renaming them with the `*.init` start controls would blur two different
  contracts: starting candidate values versus Powell direction-set geometry.

### Categorical Sample-Scaling Switch

Name:

- `scale.init.categorical.sample`

Role:

- logical switch controlling whether categorical initialization and direction
  controls are scaled by the sample-size factor documented in the bandwidth
  help pages.
- despite containing `.init`, it is a behavior switch, not a numeric starting
  value.

Recommendation:

- Keep untouched.
- Do not include in a numeric initialization rename.

Reason:

- Its risk is documentation clarity, not contract alignment with the continuous
  floor.
- Renaming it would create a broader behavior-option migration with no direct
  benefit to the scale-factor floor repair.

### `npregiv` Internal Bandwidth Initializer

Name:

- `bw.init`

Role:

- internal `npregiv` bandwidth initializer used to seed later stages from
  previously computed IV-regression bandwidths.
- appears in `R/npregiv.R` only.

Recommendation:

- Keep untouched.

Reason:

- It is IV-specific algorithm state, not a public continuous scale-factor
  search control.
- `npregiv` is lower priority than core estimators under the workspace policy
  and should not be swept into this naming tranche.

### `npRmpi` Session Initialization

Name:

- `npRmpi.init`

Role:

- exported MPI lifecycle/session entry point in `npRmpi`.
- appears throughout `npRmpi` documentation and tests.

Recommendation:

- Keep untouched.

Reason:

- It is package lifecycle API, unrelated to estimator bandwidth initialization.
- Renaming or mechanically touching it would be high-risk and out of scope.

### Test-Local State Variables

Names:

- `old.init` and similar local test variables.

Role:

- local test bookkeeping for saving/restoring options.

Recommendation:

- Ignore for naming purposes.

## Protected Do-Not-Touch List

The current continuous scale-factor naming tranche should not rename, migrate,
or change semantics for:

- `dfac.init`
- `lbd.init`
- `hbd.init`
- `initc.dir`
- `initd.dir`
- `scale.init.categorical.sample`
- `bw.init`
- `npRmpi.init`
- test-local variables such as `old.init`

## Revised Scope For The Naming Plan

The naming plan should be explicitly limited to the continuous search controls:

- `cfac.init`
- `lbc.init`
- `hbc.init`
- `scale.factor.lower.bound`

Any implementation should add a regression guard proving that the protected
names above remain present, documented where currently public, and behaviorally
unchanged after the continuous-control rename.
