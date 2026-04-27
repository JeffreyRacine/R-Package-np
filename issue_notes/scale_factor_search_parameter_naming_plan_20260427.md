# Scale-Factor Search Parameter Naming Plan

Date: 2026-04-27

## Goal

Rename the public continuous scale-factor search controls to a coherent naming
family without changing estimator semantics, defaults, optimizer behavior,
compiled objective logic, or the already repaired scale-factor floor contract.

Approved public names:

| Old name | New name | Meaning |
| --- | --- | --- |
| `cfac.init` | `scale.factor.init` | deterministic first continuous scale-factor start |
| `lbc.init` | `scale.factor.init.lower` | lower endpoint for random continuous scale-factor starts |
| `hbc.init` | `scale.factor.init.upper` | upper endpoint for random continuous scale-factor starts |
| `scale.factor.lower.bound` | `scale.factor.search.lower` | admissibility floor for continuous fixed-bandwidth search |

This is a naming/API tranche only. It must not touch manuscript files,
simulation campaign files, compiled optimizer semantics, or unrelated
initialization controls.

Companion catalogue:

```text
/Users/jracine/Development/np-master/issue_notes/init_parameter_catalogue_20260427.md
```

That catalogue is part of this plan's scope guard.

## Comfort Assessment

I am comfortable with this plan only if the implementation remains deliberately
narrow:

- strict public rename with clear old-name diagnostics;
- no estimator, optimizer, compiled-code, or default-value changes;
- canonical user-facing object metadata, with a private legacy-object reader;
- full documentation wording updated in the same tranche as the code;
- pre/post numerical parity before any positive status claim.

Without those constraints, this would be a deceptively risky rename because the
controls cross S3 dispatch, `...` forwarding, object metadata, optimizer
handoff lists, documentation, and `npRmpi` parity.

## Highest-Standards Critique Of The First Plan

The initial plan had the right conceptual names but was still too permissive
mechanically. The following issues must be corrected before implementation.

### 1. Strict Rename Versus Alias Policy Was Undecided

Leaving strict rename versus transitional aliases as a later choice invites
inconsistent behavior across selectors.

Decision:

- use a strict public rename;
- do not support old names as aliases;
- when old names are supplied, fail fast with a clear rename diagnostic rather
  than letting R produce an obscure unused-argument error or silently ignore a
  forwarded value.

Rationale:

- version `0.70-1` is intentionally a break from `0.60-*`;
- a dual-name API increases test surface and future removal risk;
- explicit diagnostics reduce collateral damage for development scripts that
  still use the old names.

### 2. Public API, Object Metadata, And Internal Optimizer Keys Were Blurred

The old plan said to preserve internal C-side option layout, which is right,
but it did not define which names should appear in public bandwidth objects.

Decision:

- public function signatures and user-facing object metadata should use the new
  canonical names;
- compiled-call option lists may retain old internal keys if that avoids C/C++
  churn;
- a narrow adapter should translate canonical R controls into existing internal
  keys at the optimizer boundary;
- object readers should tolerate pre-rename object fields only as internal
  migration support, not as public argument aliases.

Rationale:

- the user's visible API should be coherent;
- compiled logic should not be touched for a naming-only tranche;
- saved or in-session bandwidth objects should not trigger avoidable failures
  during the development transition.

### 3. `...` Forwarding Is The Main Regression Trap

Many selectors route controls through S3 methods, bootstrap constructors, and
`...`. A naïve signature rename can either drop values or produce inconsistent
error messages.

Decision:

- build a call-boundary inventory before patching;
- add one shared old-name rejection helper for public entry points and S3
  methods that receive `...`;
- add route tests proving old names error consistently at the public boundary.

Rationale:

- `...` forwarding bugs are easy to miss and hard for users to diagnose;
- a shared helper avoids selector-specific drift.

### 4. The Plan Did Not Separate Input Names From Stored Names

The floor is not just an argument; it is carried inside bandwidth objects and
used by downstream methods. Changing the stored name everywhere in one sweep is
more dangerous than changing public input names.

Decision:

- first map every read/write of the stored floor value;
- decide per field whether it is public metadata or optimizer plumbing;
- prefer one canonical stored field `scale.factor.search.lower` for user-facing
  bandwidth objects;
- allow a private compatibility reader that returns
  `scale.factor.search.lower` when present and falls back to
  `scale.factor.lower.bound` only for pre-rename objects.

Rationale:

- summaries and printed objects should match the new public API;
- internal migration support prevents avoidable breakage without keeping an old
  user-facing argument alive.
- object metadata migration must be explicitly tested because user workflows
  often reuse bandwidth objects across calls.

### 5. The Original Gate Set Was Too General

The first plan listed validation families but did not require equivalence of
raw objective traces or object structure before and after the rename.

Decision:

- run pre/post same-seed parity on representative selectors before accepting
  the rename;
- compare objective values, bandwidths, scale factors, stored floor values, and
  raw traces where available;
- require no behavioral movement except the intended names and diagnostics.

Rationale:

- a naming tranche should be almost numerically invisible.

### 6. The Protected-Name Boundary Needs Tests

The catalogue identifies protected names, but the first plan did not require a
test that they remain untouched.

Decision:

- add a protected-name guard test or static check for:
  `dfac.init`, `lbd.init`, `hbd.init`, `initc.dir`, `initd.dir`,
  `scale.init.categorical.sample`, `bw.init`, and `npRmpi.init`.

Rationale:

- this prevents the tranche from accidentally becoming a broad initializer
  rewrite.

### 7. Documentation Needs A Canonical Wording Contract

The first plan required documentation updates but did not specify what the
documentation must say. That leaves room for pages to be technically renamed
but still conceptually ambiguous.

Decision:

- every touched `.Rd` page must use the same conceptual wording for the four
  controls;
- each page must explain the difference between initialization controls and the
  search admissibility floor;
- each page must state that these are scale factors, not physical bandwidths;
- each page must state that route-specific scaling is handled internally and
  must not be reproduced by the user.

Rationale:

- the original regression came from conceptual ambiguity around floors and
  initialization;
- naming alone is not enough unless the docs make the contract unambiguous.

## Scientific And Engineering Contract

After the rename:

```r
effective.lower <- max(scale.factor.init.lower, scale.factor.search.lower)
effective.first <- max(scale.factor.init, scale.factor.search.lower)
```

and:

```r
scale.factor.init.upper >= effective.lower
```

must hold for search initialization.

The floor applies to continuous fixed-bandwidth search candidates. It is not a
floor on explicit user-supplied bandwidth objects when bandwidth search is not
being run.

No hard-coded exponent is allowed. All scale conversion must continue to use
the existing route-specific scale helpers and bandwidth-constructor machinery.

User-facing interpretation:

- users provide dimensionless continuous scale factors;
- the package converts those scale factors to route-specific bandwidths using
  the appropriate estimator/kernel/order/dimension scaling already encoded in
  the bandwidth constructors;
- `scale.factor.init` and `scale.factor.init.lower` only control starting
  values for search;
- `scale.factor.search.lower` controls admissibility of fixed-bandwidth search
  candidates;
- explicit user bandwidths remain user bandwidths when search is not run.

## Reformulated Implementation Plan

### Phase 0: Freeze And Branch

Create a feature branch from clean `np-master/master`.

Retain all diagnostics under:

```text
/Users/jracine/Development/tmp/scale_factor_search_parameter_naming_20260427
```

Do not patch `np-npRmpi` until `np-master` is committed and installed-build
validated.

### Phase 1: Boundary Inventory Before Patching

Catalogue, in a retained inventory file:

- public function signatures exposing any old continuous control name;
- S3 methods and helper constructors that receive these controls;
- `...` forwarding sites that can carry old or new names;
- bandwidth object fields that store the floor;
- internal optimizer option lists and compiled-call boundaries;
- documentation pages, examples, tests, and demos mentioning old names.

The inventory must explicitly classify each occurrence as:

- public input name;
- user-facing stored metadata;
- internal R plumbing;
- compiled optimizer boundary;
- documentation/test reference.

No source patch should begin until this classification exists.

Required inventory output:

```text
/Users/jracine/Development/tmp/scale_factor_search_parameter_naming_20260427/boundary_inventory.tsv
```

Columns:

- `repo`
- `file`
- `line`
- `old_name`
- `new_name`
- `classification`
- `planned_action`
- `risk_note`

### Phase 2: Shared Public Resolver And Old-Name Guard

Add a single shared R helper for public continuous search controls.

Responsibilities:

- accept only canonical public names;
- reject old public names with a clear rename message;
- validate scalar finite numeric controls with existing rules;
- enforce
  `scale.factor.init.upper >= max(scale.factor.init.lower,
  scale.factor.search.lower)`;
- return stable internal names for optimizer plumbing:

```r
list(
  cfac.init = ...,
  lbc.init = ...,
  hbc.init = ...,
  scale.factor.lower.bound = ...
)
```

The helper should be used at public selector boundaries, not scattered as
ad hoc per-route validation.

Old-name diagnostic wording should follow this pattern:

```text
'cfac.init' has been renamed to 'scale.factor.init'
'lbc.init' has been renamed to 'scale.factor.init.lower'
'hbc.init' has been renamed to 'scale.factor.init.upper'
'scale.factor.lower.bound' has been renamed to 'scale.factor.search.lower'
```

If multiple old names are supplied, report all of them in one error.

### Phase 3: Object-Field Migration Helper

Add a private helper for bandwidth objects:

- canonical write field: `scale.factor.search.lower`;
- legacy read fallback: `scale.factor.lower.bound`;
- conflict error if both fields are present with different values;
- no support for old argument names through this helper.

This isolates stored-object migration from public argument compatibility.

Stored metadata policy:

- new objects should write `scale.factor.search.lower`;
- summaries should display the canonical new name or a plain-language label
  tied to it;
- internal optimizer adapters may still use `scale.factor.lower.bound` as an
  implementation key;
- when reading existing bandwidth objects, the helper may use
  `scale.factor.lower.bound` only if `scale.factor.search.lower` is absent;
- if both fields are present and differ, error with a precise conflict message.

### Phase 4: `np-master` Public Surface

Update public signatures and call forwarding for core touched selectors:

- `npregbw`
- `npudensbw`
- `npudistbw`
- `npcdensbw`
- `npcdistbw`
- `npindexbw`
- `npscoefbw`
- shared/default/internal bandwidth constructors only where they are part of
  the public call path.

Hold defaults constant by route:

```r
scale.factor.init = current cfac.init default
scale.factor.init.lower = current lbc.init default
scale.factor.init.upper = current hbc.init default
scale.factor.search.lower = NULL
```

At the optimizer boundary, translate to existing internal keys if that avoids
compiled-code churn.

### Phase 5: Documentation

Update `.Rd` files in the same tranche:

- argument lists;
- Details sections;
- examples;
- error text expectations.

Each page must clearly distinguish:

- deterministic first start;
- random multistart range;
- search admissibility floor;
- explicit `bandwidth.compute = FALSE` storage exception.

The protected names from the catalogue must remain documented exactly where
they are still public.

Required documentation wording:

1. `scale.factor.init`

   - dimensionless deterministic first continuous scale-factor start for
     bandwidth search;
   - used for the first Powell/search start when continuous fixed-bandwidth
     search is run;
   - if `scale.factor.search.lower` is larger, the effective first start is
     raised to that floor;
   - not used to clamp explicit user-supplied bandwidth objects when
     `bandwidth.compute = FALSE`.

2. `scale.factor.init.lower`

   - dimensionless lower endpoint for random continuous scale-factor starts
     used in multistart bandwidth search;
   - the effective lower endpoint is
     `max(scale.factor.init.lower, scale.factor.search.lower)`;
   - not a search admissibility floor by itself.

3. `scale.factor.init.upper`

   - dimensionless upper endpoint for random continuous scale-factor starts
     used in multistart bandwidth search;
   - must be no smaller than
     `max(scale.factor.init.lower, scale.factor.search.lower)`;
   - errors if it cannot generate a feasible random-start interval.

4. `scale.factor.search.lower`

   - dimensionless lower admissibility bound for continuous fixed-bandwidth
     search candidates;
   - applies during bandwidth search and optimizer handoff/refinement;
   - raises the effective deterministic first start and random-start lower
     endpoint when needed;
   - does not override explicit user-supplied bandwidths outside search.

Required Details paragraph template:

```text
The `scale.factor.*` controls are dimensionless search controls. The package
converts scale factors to bandwidths using the estimator-specific scaling
already encoded in the bandwidth object, including the kernel order and the
number of continuous variables relevant for that estimator. Users should not
pre-multiply these controls by sample-size or standard-deviation factors.
`scale.factor.init` controls the deterministic first search start,
`scale.factor.init.lower` and `scale.factor.init.upper` control the random
multistart interval, and `scale.factor.search.lower` is the lower admissibility
bound for continuous fixed-bandwidth search candidates.
```

Each touched help page must also include a sentence distinguishing
`scale.factor.search.lower` from categorical initialization controls:

```text
Categorical search-start controls such as `dfac.init`, `lbd.init`, and
`hbd.init` have separate semantics and are not affected by
`scale.factor.search.lower`.
```

### Phase 6: `np-master` Tests

Add or update focused tests proving:

- canonical names are accepted on every touched selector;
- old names error with clear rename diagnostics;
- conflicting object fields error;
- legacy object field fallback works only for object metadata migration;
- bad `scale.factor.init.upper` errors against the effective lower endpoint;
- explicit bandwidth objects remain unclamped when search is not run;
- final fixed candidates still obey `scale.factor.search.lower`;
- protected names remain present and unchanged;
- no hard-coded exponent is introduced.

Add installed-build sentinels covering at least:

- one simple route;
- one representative heavier `npcdens` route;
- `npindexbw`;
- `npscoefbw`.

Documentation tests/checks:

- every touched `.Rd` file contains all four canonical names;
- no touched `.Rd` file lists old continuous names as public arguments;
- protected categorical and direction-set names remain where currently public;
- examples use canonical names only.

### Phase 7: Numerical And Structural Parity

Before accepting the rename, compare pre/post installed builds using identical
seeds and data.

Required comparisons:

- objective values;
- selected bandwidths;
- selected scale factors;
- stored floor metadata;
- raw objective traces where available;
- printed summaries for the expected name changes.

Acceptance rule:

- numerical behavior must be unchanged within existing route tolerances;
- differences are allowed only for public names, stored metadata names,
  diagnostic messages for old arguments, and documentation text.

Stop rule:

- if any route changes selected bandwidths, scale factors, objective values, or
  raw trace values beyond established tolerance, stop and diagnose before
  proceeding to documentation polish or `npRmpi`.
- if the only difference is object field naming, verify that downstream
  estimator, `summary`, `predict`, and adjacent methods still consume the
  object correctly.

### Phase 8: `np-master` Installed Gates And Commit

Run:

- focused naming tests;
- existing scale-factor floor tests;
- protected-name static/test guard;
- documentation-name checks;
- sentinel scripts from Phase 6;
- installed tarball checks from retained tarball;
- `R CMD check --as-cran`.

Commit only after these gates pass.

Suggested commit shape:

1. resolver, object-field migration, and R surface;
2. docs, tests, and retained validation notes if the first commit is too large.

### Phase 9: Surgical `npRmpi` Port

Only after `np-master` is committed and green:

- create a separate `npRmpi` branch;
- port by surgical patch, never file copy;
- preserve MPI lifecycle and helper policy;
- preserve runtime independence; no `np::` bridge calls;
- keep `npRmpi.init` untouched.

Run `npRmpi` gates:

- focused naming tests through shared MPI test helpers;
- existing scale-factor floor tests;
- sentinel matrix with `nslaves = 1` and `nslaves = 2`;
- installed tarball checks;
- `R CMD check --as-cran`.

Commit only after gates pass.

### Phase 10: Merge Discipline

Merge only after both repos are independently green:

- feature branch into `np-master/master`;
- feature branch into `np-npRmpi/npRmpi`.

Do not push unless explicitly requested.

## Explicit Non-Goals

Do not rename or change:

- `dfac.init`
- `lbd.init`
- `hbd.init`
- `initc.dir`
- `initd.dir`
- `scale.init.categorical.sample`
- `bw.init`
- `npRmpi.init`

Do not change:

- defaults;
- estimator semantics;
- compiled objective functions;
- optimizer handoff behavior;
- simulation harnesses or campaign files;
- manuscript files.
