# Cross-Package `.Rd` User-Model Documentation Plan

Date: 2026-04-29

Scope:
- `/Users/jracine/Development/np-master`
- `/Users/jracine/Development/np-npRmpi`
- `/Users/jracine/Development/crs`
- `/Users/jracine/Development/spcovar`

Status: plan only. Do not edit package code, `.Rd`, tests, or generated
artifacts until this plan is accepted for execution.

Related documents:
- `/Users/jracine/Development/np-master/issue_notes/scale_factor_search_parameter_naming_plan_20260427.md`
- `/Users/jracine/Development/np-master/issue_notes/rd_documentation_style_ledger_20260428.md`
- `/Users/jracine/Development/np-master/issue_notes/release_readiness_arg_rd_cleanup_plan_20260428.md`

## Goal

Improve the user-facing `.Rd` documentation across `np`, `npRmpi`, `crs`, and
`spcovar` by organizing large argument blocks around the user's mental model
rather than historical accretion order. Preserve the current successful
subsection-header style and refine it where necessary; do not discard it for any
function.

Within each mental-model group, order arguments alphabetically by argument name.
If alphabetical order feels wrong because the group contains controls with
different user purposes, treat that as evidence that the group is too broad and
split it into clearer mental-model groups. Reserve non-alphabetical ordering for
rare, explicitly documented exceptions where a very small sequence is itself
the concept being taught.

For `np` and `npRmpi`, this includes a special focus on the continuous
`scale.factor.*` search controls. For `crs` and `spcovar`, apply the same
mental-model grouping principle only where the page has enough arguments to
benefit; do not invent scale-factor wording where those controls are not public.

This tranche is documentation-only unless a later explicit decision changes the
scope. The current code contract should be preserved:

- `scale.factor.init = 0.5`
- `scale.factor.init.lower = 0.1`
- `scale.factor.init.upper = 2.0`
- `scale.factor.search.lower = NULL`

The `NULL` default for `scale.factor.search.lower` is intentional. It means:

1. inherit a stored value from an existing bandwidth object when one is
   available;
2. otherwise resolve to the package default search floor, currently `0.1`.

Do not change this formal default merely for cosmetic consistency.

## Protected Documentation Style

The existing subsection-header approach is now a protected readability feature
for this documentation campaign.

Rules:

1. Preserve useful existing `\subsection{...}{...}` groupings.
2. Refine overbroad headings into clearer mental-model headings only when the
   rendered help becomes easier to scan.
3. Do not collapse grouped pages back into flat argument lists.
4. Do not remove subsection headers merely to simplify editing.
5. If a page already has excellent subsection organization, keep the structure
   and make only local ordering/wording improvements.
6. Every `\subsection{...}{...}` inside `\arguments{}` must keep nonempty body
   text and all real arguments must remain top-level `\item{arg}{...}` siblings.

## Documentation Problem

The current grouped `.Rd` structure is much better than the historical flat
argument lists, but large sections such as `Search Initialization, Kernels, And
Support` can still interleave unrelated controls:

- continuous scale-factor starts,
- continuous direction-set controls,
- categorical starts and direction-set controls,
- kernel choices,
- support bounds,
- tolerances,
- penalties.

This makes the `scale.factor.*` contract hard to read because the four controls
are split apart:

- `scale.factor.init`
- `scale.factor.init.lower`
- `scale.factor.init.upper`
- `scale.factor.search.lower`

These should be documented contiguously and with defaults.

The same general problem can occur outside `np`/`npRmpi`: options added over
time can become locally disordered even when the page now has good subsection
headings. The cross-package correction is not to rewrite every page, but to
audit rendered help pages and improve only those where user comprehension is
materially improved by:

- preserving high-level positional/main inputs at the top;
- grouping optional controls by mental model;
- ordering controls alphabetically within each group;
- splitting any group whose contents are too broad to scan alphabetically;
- documenting defaults consistently;
- explaining formal `NULL` defaults by also naming the effective package
  default or inheritance rule when the code supplies one;
- avoiding runtime/API changes.

## User-Model Ordering Contract

Within large help pages, organize optional controls by conceptual clusters
first. Within each conceptual cluster, use alphabetical order by argument name.
If an argument group appears to need a non-alphabetical internal order to be
understandable, first split it into more precise groups.

General ordering rule:

1. Keep user entry-point arguments first under a combined leading subsection
   when present:
   `Data, Bandwidth Inputs And Formula Interface`.
   This section collects immutable/main call-pattern arguments users plausibly
   pass positionally:
   - `formula`, `bws`, `data`;
   - `txdat`, `tydat`, `tzdat`, `tdat`;
   - `exdat`, `eydat`, `ezdat`, `edat`;
   - package-specific object or model inputs.
2. Keep `Additional Arguments` as the final subsection when present.
3. Between the leading entry-point subsection and final `Additional Arguments`,
   order subsection labels alphabetically unless a package-specific user
   workflow clearly argues otherwise and the exception is recorded.
4. Within each subsection, order arguments alphabetically by argument name.
5. Treat proposed non-alphabetical ordering as a warning sign. Prefer splitting
   an overbroad subsection into smaller user-model groups rather than relying on
   hidden editorial ordering.
6. Use a `Miscellaneous Controls` or `Additional Controls` subsection only as a
   last resort for a small number of genuinely unrelated arguments; if the
   miscellaneous group grows large, split it again.
7. Allow a non-alphabetical exception only when a very small sequence is itself
   a public contract or workflow that would be harder to understand if sorted.
   Record the exception in the page ledger.

## Default-Value Documentation Contract

Argument documentation should make defaults easy for users to discover without
forcing them to inspect R source.

Rules:

1. For every touched argument with a formal default in `\usage{}`, state the
   default in the argument text unless the default is already clear in an
   immediately adjacent table, list, or parent sentence.
2. If the formal default is `NULL` but the package resolves `NULL` through an
   option list, bandwidth object, helper default, estimator family, or other
   documented internal default, state both facts:
   - formal default: `NULL`;
   - effective behavior: inherit or use the named package default when
     applicable.
3. If `NULL` means "not supplied", "infer from data", or "skip this feature",
   say that plainly rather than inventing a numeric or character default.
4. If the default depends on another argument, data dimension, method family, or
   package option, document the rule rather than freezing one misleading value.
   Examples include `min(2,ncol(xdat))`, method-specific fallbacks, and
   option-backed controls.
5. Verify defaults against live formals, option objects, or the relevant helper
   resolver before editing prose. Do not infer defaults solely from old `.Rd`
   wording.
6. Do not change R formals, option defaults, runtime behavior, object storage,
   or examples to make the documentation look more uniform.
7. Record any default whose source is non-obvious in the page ledger so later
   reviewers can verify it without repeating the excavation.

For `np`/`npRmpi` bandwidth-selector pages, the four continuous scale-factor
controls are already readable in alphabetical order, which also preserves the
natural user story:

1. `scale.factor.init`
2. `scale.factor.init.lower`
3. `scale.factor.init.upper`
4. `scale.factor.search.lower`

Preferred order for relevant `np`/`npRmpi` bandwidth-selector groups:

1. Continuous scale-factor search starts
   - `scale.factor.init`
   - `scale.factor.init.lower`
   - `scale.factor.init.upper`
   - `scale.factor.search.lower`
2. Continuous direction-set controls
   - `cfac.dir`
   - `dfc.dir`
   - `hbc.dir`
   - `initc.dir`
   - `lbc.dir`
3. Categorical search starts and direction-set controls
   - `dfac.dir`
   - `dfac.init`
   - `hbd.dir`
   - `hbd.init`
   - `initd.dir`
   - `lbd.dir`
   - `lbd.init`
   - `scale.init.categorical.sample`
4. Kernel type controls
   - `ckerorder`
   - `ckertype`
   - `okertype`
   - `ukertype`
   - family-specific analogues such as `cxkertype`, `cykertype`, `oxkertype`,
     `oykertype`, `uxkertype`, and `uykertype` where present
5. Continuous kernel support controls
   - `ckerbound`
   - `ckerlb`
   - `ckerub`
   - family-specific bounded-support analogues such as `cxkerbound`,
     `cxkerlb`, `cxkerub`, `cykerbound`, `cykerlb`, and `cykerub` where present
6. Numerical search and tolerance controls
   - `ftol`
   - `invalid.penalty`
   - `itmax`
   - `nmulti`
   - `penalty.multiplier`
   - `remin`
   - `small`
   - `tol`
   - `transform.bounds`
   - family-specific search controls that naturally belong with numerical
     tolerances

This ordering should be adapted by page. Do not force controls into a page where
they are not already public for that function family.

For `crs` and `spcovar`, derive analogous groups from each page's actual
arguments. Examples may include:

- data/model inputs;
- bandwidth or smoothing controls;
- kernel/basis controls;
- optimization and numerical tolerance controls;
- prediction/evaluation controls;
- plotting/interval controls;
- compatibility and additional arguments.

The rule is not to force `np` headings into non-`np` packages; the rule is to
preserve useful subsection headers and make the user's workflow easy to scan.

## Example Replacement For The Broad Search/Kernels Bucket

When a page currently has a broad group such as:

```rd
\subsection{Search Initialization, Kernels, And Support}{
  These controls set numerical search starts, kernel choices, and support
  bounds.
}
```

replace it on large pages with sharper user-model headings:

```rd
\subsection{Continuous Scale-Factor Search Initialization}{
  These controls define deterministic and random continuous scale-factor starts
  and the lower admissibility floor for fixed-bandwidth search.
}

\subsection{Continuous Direction-Set Search Controls}{
  These controls set Powell direction-set initialization for continuous
  variables.
}

\subsection{Categorical Search Initialization}{
  These controls set categorical search starts and categorical direction-set
  initialization.
}

\subsection{Kernel Type Controls}{
  These controls choose continuous, unordered, and ordered kernels.
}

\subsection{Continuous Kernel Support Controls}{
  These controls choose and parameterize bounded support for continuous kernels.
}

\subsection{Numerical Search And Tolerance Controls}{
  These controls set optimizer tolerances, restart behavior,
  invalid-candidate penalties, and bounded search transformations.
}
```

Every `\subsection{...}{...}` inside `\arguments{}` must have nonempty body
text, and all real arguments must remain top-level `\item{arg}{...}` siblings
after the relevant heading. Do not nest `\item{}` entries inside the subsection
body.

## Preferred `.Rd` Heading Strategy

Use real empty sibling `\subsection{...}{}` headings inside `\arguments{}` where
the page is large enough to benefit.

For large pages, split the existing broad section into clearer headings when it
improves scanability:

- `Continuous Scale-Factor Search Initialization`
- `Continuous Direction-Set Search Controls`
- `Categorical Search Initialization`
- `Kernel Type Controls`
- `Continuous Kernel Support Controls`
- `Numerical Search And Tolerance Controls`

For moderate pages, it is acceptable to keep a broader heading such as
`Search Initialization, Kernels, And Support` if the internal order follows the
user-model contract above and the page remains readable.

Do not use fake group headings as `\item{...}{...}`.
Do not use a broad heading on pages where it leaves more than about ten
heterogeneous arguments in one block; split those pages into the sharper
subsections above.

## Required Scale-Factor Wording

Use consistent wording across `np` and `npRmpi`, adapted only where a function
family requires estimator-specific wording.

### `scale.factor.init`

Suggested wording:

> Deterministic initial scale factor for continuous fixed-bandwidth search.
> Defaults to `0.5`. The value supplied by the user is not rewritten, but the
> effective first start passed to the optimizer is
> `max(scale.factor.init, scale.factor.search.lower)`. See Details.

### `scale.factor.init.lower`

Suggested wording:

> Lower endpoint for random continuous scale-factor starts. Defaults to `0.1`.
> The value supplied by the user is not rewritten, but the effective random-start
> lower endpoint is `max(scale.factor.init.lower,
> scale.factor.search.lower)`. See Details.

### `scale.factor.init.upper`

Suggested wording:

> Upper endpoint for random continuous scale-factor starts. Defaults to `2.0`.
> It must be greater than or equal to the effective lower endpoint,
> `max(scale.factor.init.lower, scale.factor.search.lower)`; otherwise bandwidth
> search errors rather than silently expanding the interval. See Details.

### `scale.factor.search.lower`

Suggested wording:

> Optional nonnegative scalar giving the hard lower admissibility bound for
> continuous fixed-bandwidth search candidates. Defaults to `NULL`. If `NULL`,
> an existing bandwidth object's stored value is inherited when available;
> otherwise the package default `0.1` is used. This floor applies to
> computed/search bandwidth candidates and to effective search starts only. It
> does not rewrite explicit bandwidths supplied for storage with
> `bandwidth.compute = FALSE`. Final fixed-bandwidth search candidates must also
> have a finite valid raw objective value.

## Required Details Paragraph

Each touched bandwidth-selector page that exposes the `scale.factor.*` controls
should contain a Details paragraph with this substance:

```text
The scale.factor.* controls are dimensionless search controls. The package
converts scale factors to bandwidths using the estimator-specific scaling
encoded in the bandwidth object, including kernel order and the number of
continuous variables relevant for the estimator. Users should not pre-multiply
these controls by sample-size or standard-deviation factors.

scale.factor.init controls the deterministic first search start.
scale.factor.init.lower and scale.factor.init.upper define the random
multistart interval. scale.factor.search.lower is the lower admissibility bound
for continuous fixed-bandwidth search candidates. The effective first start is
max(scale.factor.init, scale.factor.search.lower), and the effective
random-start lower endpoint is max(scale.factor.init.lower,
scale.factor.search.lower). scale.factor.init.upper must be at least that
effective lower endpoint; the package errors rather than silently expanding the
user's interval.

When scale.factor.search.lower is NULL, an existing bandwidth object's stored
floor is inherited when available; otherwise the package default 0.1 is used.
Explicit bandwidths supplied for storage with bandwidth.compute = FALSE are not
rewritten by the search floor.

Categorical search-start controls such as dfac.init, lbd.init, and hbd.init
have separate semantics and are not affected by scale.factor.search.lower.
```

Convert to valid `.Rd` markup with `\code{}` and `\eqn{}` as appropriate.

## Candidate Pages

Audit all four repos for pages that would materially benefit from mental-model
argument grouping or ordering.

For `np` and `npRmpi`, audit every page exposing any `scale.factor.*` argument.

Expected high-priority pages include:

- `man/np.regression.bw.Rd`
- `man/np.density.bw.Rd`
- `man/np.distribution.bw.Rd`
- `man/np.condensity.bw.Rd`
- `man/np.condistribution.bw.Rd`
- `man/np.smoothcoef.bw.Rd`
- `man/np.singleindex.bw.Rd`
- `man/np.plregression.bw.Rd`

Apply the same audit to the corresponding `npRmpi` pages, preserving package
specifics.

For `crs` and `spcovar`, first inventory pages with:

- long `\arguments{}` sections;
- existing broad/mixed subsection headings;
- option blocks that have clearly grown organically;
- user-facing controls whose defaults or semantics are not easy to scan.

Do not edit short, already-clear pages merely for uniformity.

## Initial Execution Outline Superseded By Reformulated Plan

This outline records the first-pass idea for the `np`/`npRmpi` scale-factor
focus. Use the `Reformulated Plan` below for execution across all four
packages; it adds the risk controls, pilot gates, rendered-help review, and
cross-repo discipline required for release-hardening work.

1. Inventory
   - List every `np`/`npRmpi` `.Rd` page containing `scale.factor.init`,
     `scale.factor.init.lower`, `scale.factor.init.upper`, or
     `scale.factor.search.lower`.
   - For each page, record whether it already uses `\arguments{}` subheaders.
   - Record whether the page has all four controls or only a subset.

2. Pilot one page
   - Start with `np-master/man/np.regression.bw.Rd`.
   - Replace the broad `Search Initialization, Kernels, And Support` group with
     the sharper subsections where appropriate.
   - Reorder only the relevant `.Rd` argument entries by user-model group and
     alphabetical order within group.
   - Do not touch R function formals or `\usage{}`.
   - Add defaults and the explicit `NULL` inheritance wording.
   - Render/check the page with `tools::checkRd()`.

3. Review rendered help
   - Inspect the rendered help output, not only source text.
   - Confirm the user can quickly answer:
     - what the deterministic start is;
     - what the random start interval is;
     - what the search floor is;
     - why `scale.factor.search.lower = NULL` is intentional;
     - what happens when the floor exceeds `scale.factor.init.lower`;
     - what happens when the floor exceeds `scale.factor.init.upper`.

4. Expand within `np-master`
   - Apply the same pattern family by family.
   - Keep edits mechanical and documentation-only.
   - On large pages, split overbroad mixed-purpose groups into the user-model
     subsections above.
   - Do not reorder unrelated argument groups unless they are part of the same
     reader-model cleanup.

5. Port surgically to `npRmpi`
   - Do not copy files between repos.
   - Apply equivalent `.Rd` edits by page, preserving `npRmpi`-specific wording,
     MPI examples, and package names.

6. Validation
   - Run `tools::checkRd()` on every touched `.Rd` file.
   - Run `R CMD build` for each package after `.Rd` edits.
   - Run targeted `R CMD check` documentation stages or full check if the build
     reveals documentation warnings.
   - Grep gates:
     - no touched page says the user value itself is rewritten;
     - no touched page says `scale.factor.search.lower` has formal default
       `0.1`;
     - all touched pages that expose all four controls list them contiguously;
     - all touched large pages use mental-model subsections instead of one
       overbroad search/kernels/support bucket;
     - arguments are alphabetical within each subsection unless a rare,
       ledger-recorded exception is justified after attempting to split the
       group;
     - all touched pages mention defaults.

## Engineering Critique Of This Plan

The current plan has the right user-facing goal, but several risks need to be
controlled before execution.

1. `.Rd` reordering can drift away from function behavior.
   - The plan correctly says not to change R formals or `\usage{}`, but it
     should also require an explicit unchanged-usage gate. Otherwise a future
     edit could make the help look cleaner while accidentally changing the
     visible call signature.

2. Strict alphabetical ordering inside a well-scoped group is helpful for
   lookup, but it can expose that a group is actually too broad.
   - If an editor wants to reorder a group semantically, the safer user-facing
     response is to split the group into clearer subsections first.
   - Lower/upper pairs can be kept together by choosing better subsection
     boundaries and wording; do not use hidden non-alphabetical ordering as the
     primary organizing principle.
   - The four `scale.factor.*` controls work well in alphabetical order and
     should remain contiguous.

3. The plan could over-split moderate pages.
   - Too many `\subsection{}` headings can make smaller help pages feel choppy.
     The split should be mandatory only for large mixed blocks; shorter pages
     may keep broader headings if the internal order is clear.

4. The plan does not yet require a rendered-output review artifact.
   - `tools::checkRd()` proves syntax, not readability. The campaign goal is
     user comprehension, so each pilot should retain rendered help text or a
     small before/after excerpt under `/Users/jracine/Development/tmp`.

5. The plan risks doc-only false confidence.
   - Even documentation-only edits can break package checks through malformed
     Rd markup, stale aliases, or examples. Build/check gates must be run after
     each repo tranche, not only at the end.

6. Cross-repo parity needs a diff discipline.
   - `np` and `npRmpi` must be semantically aligned, but direct file copying is
     prohibited and unsafe. The plan should require a page-by-page parity
     checklist rather than textual mirroring.

7. Cross-package generalization can dilute the core win.
   - `crs` and `spcovar` should receive the same readability treatment where
     needed, but their argument taxonomies differ. Forcing `np` headings onto
     them would create a new kind of documentation noise.

8. Scope creep is likely.
   - This tranche should not become a general documentation rewrite. The target
     is argument organization, subsection preservation/refinement, and
     scale-factor wording where those controls are public. Other doc problems
     should be logged separately unless they block this readability improvement.

9. The new default-documentation rule could become too broad.
   - Requiring a full source excavation for every argument in every touched page
     would slow the campaign, invite mistakes, and increase collateral edits.
   - Refinement: treat default documentation as a touched-argument gate, not a
     mandate to rewrite every argument. Verify and state defaults for arguments
     whose entries are moved, regrouped, renamed by heading, or whose existing
     prose is visibly unclear. Log broader default gaps separately unless they
     directly affect the touched user-model groups.

10. Combining formula/data inputs with default-method data inputs can blur
    interface roles.
    - A single leading section is good for users, but it must not imply that
      every argument applies to every S3 method.
    - Refinement: the leading section wording should say these are the primary
      data, bandwidth, and formula entry points. Individual entries should
      preserve method-specific meaning where needed, and `\usage{}` remains the
      source of truth for exact method signatures.

11. Alphabetizing middle subsection labels can separate workflow-neighbor
    concepts.
    - This improves scanability, but it could place a setup concept far away
      from a dependent control.
    - Refinement: use cross-references in prose (`See ...`) rather than
      violating alphabetical order. If the split makes rendered help worse,
      record the page as deferred instead of forcing the global rubric.

## Final Engineering Critique And Refinements

This final pass adds the remaining controls needed to maximize the chance of a
clean documentation win with minimal regression risk.

1. The plan must protect against accidental argument loss or duplication.
   - Reordering `.Rd` entries by hand can drop an argument, duplicate an entry,
     or move an argument away from the only page that explains it.
   - Refinement: every touched page needs a before/after argument-set invariant:
     same argument names, same aliases, same `\usage{}` block, unless an
     explicit page note explains a deliberate documentation-only correction.

2. The plan should not rely only on manual pre-sweep judgment.
   - Deep understanding is required, but basic invariants should be generated
     mechanically so the human review can focus on semantics.
   - Refinement: create small extraction helpers under the scratch root to
     produce per-page argument lists, subsection headings, and usage snapshots
     before and after edits.

3. The plan still risks touching live release trees too early.
   - Documentation-only work can still break builds. Given release proximity,
     the pilot should earn trust before edits are applied to live repos.
   - Refinement: perform pilot and broad inventories in a detached worktree or
     scratch copy first when feasible; only apply to live repos after the pilot
     page has passed parse, rendered-help, invariant, and build checks.

4. Rendered help must be the review surface.
   - Source ordering can look clean while rendered help is too choppy or
     verbose.
   - Refinement: for each pilot/batch, retain `Rd2txt` rendered output before
     and after. For at least one representative large page per package, inspect
     rendered help manually before continuing.

5. Package-specific taxonomies must remain package-specific.
   - The cross-package scope is helpful, but `crs` and `spcovar` should not be
     forced into the `np` bandwidth-selector mental model.
   - Refinement: `crs` and `spcovar` require their own classification ledgers
     and package-specific group names. Reuse style principles, not headings.

6. Good existing subsection organization should be treated as a keep unless
   evidence says otherwise.
   - The new subsection approach is a readability improvement and should not be
     churned unnecessarily.
   - Refinement: page notes must explicitly classify each existing subsection
     as `keep`, `split`, `rename`, or `skip`; default is `keep`.

7. The tranche needs stopping rules.
   - If the first pilot makes the page more complex or brittle, expanding would
     multiply risk.
   - Refinement: after each package pilot, stop if rendered help is not clearly
     easier to scan, if invariants fail, or if check/build warnings appear.
     Repair or defer rather than pushing onward.

## Reformulated Plan

Use this reformulated plan for execution. It supersedes the earlier
high-level execution outline while preserving its intent.

### Phase 0: Safety Baseline

1. Confirm all four repos' status before edits.
2. Record the exact list of candidate `.Rd` pages in all four repos.
3. Save pre-edit snapshots of:
   - `\usage{}` blocks for candidate pages;
   - the relevant `\arguments{}` blocks;
   - `tools::Rd_db()` parse success for candidate pages.
4. Create a scratch root under:
   `/Users/jracine/Development/tmp/cross_package_rd_user_model_20260429`
5. Do not edit R code, C code, tests, examples, vignettes, or `\usage{}` blocks
   in this tranche.
6. Record existing subsection headings for every candidate page so useful
   headings are preserved unless the rendered help clearly improves by refining
   them.
7. For pilot work, prefer a detached worktree or scratch copy until the pilot
   page passes parse, rendered-help, invariant, and build checks. Apply to the
   live repo only after the pilot earns keep status.
8. Create extraction helpers under the scratch root to capture:
   - `.Rd` argument names in order;
   - subsection headings in order;
   - `\usage{}` blocks;
   - rendered `Rd2txt` help output.

### Phase 1: Pre-Sweep And Argument Classification

Do not begin `.Rd` edits until this pre-sweep is complete for the package being
edited.

For each package, create a classification ledger under the scratch root before
touching any `.Rd` file:

```text
/Users/jracine/Development/tmp/cross_package_rd_user_model_20260429/<package>_argument_classification.tsv
```

The ledger should record, at minimum:

- `.Rd` file;
- documented function or alias;
- current subsection heading;
- argument name;
- whether the argument appears in `\usage{}`;
- whether it appears in the corresponding R function formals where applicable;
- whether it is a likely positional/main input;
- proposed mental-model group;
- proposed within-group order;
- whether ordering is alphabetical or semantic-exception;
- subsection action: `keep`, `split`, `rename`, `new`, or `skip`;
- notes on unclear purpose or package-specific behavior.

Pre-sweep method:

1. Extract all candidate `.Rd` `\arguments{}` entries.
2. Extract `\usage{}` blocks and corresponding R formals where feasible.
3. Preserve leading/main positional arguments as protected:
   - formulas, bandwidth objects, data/model objects, training/evaluation data,
     and package-specific object inputs.
4. Classify optional controls by function-specific meaning, not by name alone.
5. For unclear arguments or non-obvious defaults, inspect the R implementation,
   package option objects, and relevant resolver helpers before assigning a
   group or adding default prose.
6. If an argument's purpose remains unclear after source inspection, leave it in
   its current relative location and note it in the ledger rather than guessing.
7. Use the ledger to decide whether a page should be edited, skipped, or
   deferred.
8. Default to preserving existing subsection headings unless the ledger records
   why a split or rename improves rendered help.

This phase is intentionally function-by-function. The campaign should not rely
on global text replacement or a purely alphabetical reorder.

### Phase 2: Pilot `np.regression.bw.Rd`

Pilot exactly one page:

- `/Users/jracine/Development/np-master/man/np.regression.bw.Rd`

Pilot rules:

1. Preserve the existing `\usage{}` block byte-for-byte unless a parse defect is
   discovered and separately approved.
2. Reorganize only the affected `\arguments{}` entries.
3. Replace the broad `Search Initialization, Kernels, And Support` subsection
   with sharper mental-model subsections only where it improves rendered
   scanability.
4. Use nonempty `\subsection{Title}{one sentence}` bodies and keep real
   arguments as top-level `\item{arg}{...}` siblings.
5. Add defaults directly to the four `scale.factor.*` entries.
6. Explain `scale.factor.search.lower = NULL` as:
   - inherit a stored bandwidth-object floor when available;
   - otherwise resolve to package default `0.1`.
7. Do not claim the user-supplied values are rewritten. Say only that effective
   optimizer starts/endpoints are computed from the user values and the search
   floor.
8. For every other touched argument whose default is visible in `\usage{}` or
   supplied by the relevant package option/resolver, add or preserve clear
   default wording when the rendered help would otherwise require source
   inspection.

Pilot ordering rules:

1. `Data, Bandwidth Inputs And Formula Interface` comes first when those
   arguments are present.
2. `Additional Arguments` comes last when present.
3. Middle subsection labels are alphabetized by heading.
4. Within groups, order alphabetically by argument name.
5. If alphabetical order makes a group feel awkward, split the group before
   considering a non-alphabetical exception.
6. Use a miscellaneous/additional group only for a small tail of genuinely
   unrelated controls; large miscellaneous groups fail the readability goal.
7. Record any remaining non-alphabetical exception in the page audit notes, with
   the reason and the exact argument sequence.

Pilot validation:

1. `tools::checkRd()` on the touched page.
2. Render the help text and save a before/after excerpt under the scratch root.
3. Grep gates on the touched page:
   - no legacy `cfac.init`, `lbc.init`, `hbc.init`, or
     `scale.factor.lower.bound` listed as public arguments;
   - no phrase implying the user value itself is rewritten;
   - `scale.factor.search.lower` does not claim formal default `0.1`;
   - all four `scale.factor.*` controls are contiguous.
4. Confirm `\usage{}` is unchanged.
5. Confirm before/after `.Rd` argument names are identical unless an explicit
   page note documents a deliberate correction.
6. Confirm no argument name appears twice after the edit.
7. Confirm default wording for touched arguments matches the live formal,
   package option, or documented resolver behavior.
8. Confirm the subsection action ledger classifies each existing heading as
   `keep`, `split`, `rename`, or `skip`.

Stop after the pilot if rendered help is not clearly better.

### Phase 3: Expand Within `np-master`

Apply the accepted pilot pattern page by page, not by bulk replacement.

For each candidate page:

1. Determine whether it is large enough to split into sharper subsections.
2. Preserve `\usage{}` and R function formals.
3. Reorder only `\arguments{}` entries that belong to the target
   search/kernel/support blocks.
4. Keep page-specific controls in the nearest appropriate mental-model group.
5. Avoid moving unrelated high-level inputs such as `formula`, `bws`, `data`,
   `txdat`, `tydat`, `tdat`, `edat`, `exdat`, or `eydat` away from their
   current main-input positions.
6. Run `tools::checkRd()` immediately after each page or small family batch.
7. Append page notes to the scratch audit ledger:
   - groups introduced;
   - non-alphabetical semantic exceptions;
   - scale-factor wording status;
   - subsection actions;
   - before/after argument-set invariant result;
   - rendered-help artifact path;
   - validation command and result.

### Phase 4: Build Gate For `np-master`

After all selected `np-master` pages are edited:

1. Run the grep gates across touched pages.
2. Run before/after invariant checks across touched pages:
   - same `\usage{}` blocks;
   - same `.Rd` argument names unless explicitly documented;
   - no duplicated argument entries.
3. Run `tools::checkRd()` across touched pages.
4. Run `R CMD build`.
5. If build succeeds, run at least the documentation-relevant portion of
   `R CMD check`; if any warning/note appears attributable to this tranche,
   fix before touching `npRmpi`.

Do not port to `npRmpi` until `np-master` passes these gates.

### Phase 5: Surgical `npRmpi` Port

Port the accepted documentation pattern page by page.

Rules:

1. Do not copy `.Rd` files between repos.
2. Preserve `npRmpi`-specific text, examples, aliases, MPI references, and
   package names.
3. Match the semantic organization and scale-factor wording from `np-master`
   where the public controls are equivalent.
4. Record any deliberate divergence in the scratch audit ledger.
5. Run `tools::checkRd()` after each page or small family batch.

### Phase 6: Build Gate For `npRmpi`

After all selected `npRmpi` pages are edited:

1. Run the same grep gates across touched pages.
2. Run `tools::checkRd()` across touched pages.
3. Run `R CMD build`.
4. Run the documentation-relevant `R CMD check` gate, escalating to full check
   if any documentation warning, namespace warning, or example warning appears.

### Phase 7: `crs` Documentation Pass

After `np` and `npRmpi` pass their gates, perform a separate `crs` pass.

Rules:

1. Start with inventory only; identify pages where argument ordering or broad
   option blocks materially impair readability.
2. Preserve existing useful subsection headers.
3. Refine headings only when the rendered help improves.
4. Use the leading `Data, Bandwidth Inputs And Formula Interface` section where
   applicable, keep `Additional Arguments` last when present, and preserve
   `\usage{}`.
5. Alphabetize middle subsection labels and argument names within
   mental-model groups.
6. If a group seems to need semantic ordering, split or defer rather than
   silently using non-alphabetical order.
7. Do not import `np`-specific headings or scale-factor wording unless the page
   actually exposes analogous controls.
8. Use a `crs`-specific classification ledger and group names.
9. Run `tools::checkRd()` after each page or small batch, then `R CMD build`
   and the documentation-relevant `R CMD check` gate for `crs`.

### Phase 8: `spcovar` Documentation Pass

After `crs` passes its gates, perform a separate `spcovar` pass.

Rules:

1. Because `spcovar` is expected to be simpler, default to minimal edits.
2. Preserve any existing good subsection-header organization.
3. Improve only pages where rendered help is clearly more readable after
   mental-model grouping or ordering.
4. Use the leading `Data, Bandwidth Inputs And Formula Interface` section where
   applicable, keep `Additional Arguments` last when present, and preserve
   `\usage{}`.
5. Alphabetize middle subsection labels and argument names within
   mental-model groups.
6. Use a `spcovar`-specific classification ledger and group names.
7. Run `tools::checkRd()` after each page or small batch, then `R CMD build`
   and the documentation-relevant `R CMD check` gate for `spcovar`.

### Phase 9: Final Cross-Package Closeout

Before declaring success:

1. Confirm every touched package has a retained pre-sweep classification ledger.
2. Compare candidate page inventories across all four repos.
3. Confirm `np`/`npRmpi` parity for equivalent pages and record any deliberate
   divergence.
4. Confirm every `np`/`npRmpi` page exposing the four `scale.factor.*` controls:
   - lists defaults;
   - explains `NULL` inheritance/default behavior;
   - keeps the four controls contiguous;
   - uses mental-model grouping appropriate to page size.
5. Confirm `crs` and `spcovar` pages preserve useful subsection headers and
   improve only pages where rendered help readability benefits.
6. Confirm no runtime files changed.
7. Confirm every touched page has:
   - unchanged `\usage{}` block;
   - unchanged argument-name set unless explicitly documented;
   - no duplicate argument entries;
   - retained rendered before/after help output.
8. Retain the scratch audit ledger and validation logs under
   `/Users/jracine/Development/tmp/cross_package_rd_user_model_20260429`.
9. Summarize any deferred documentation issues separately; do not fix them
   inside this tranche unless they directly block the scale-factor/user-model
   readability goal.

## Non-Goals

- Do not change runtime defaults.
- Do not change R function formal ordering.
- Do not change `\usage{}` ordering for visual tidiness.
- Do not rename any additional arguments.
- Do not add Gallery URLs or any external links.
- Do not remove useful existing subsection headers.
- Do not force `np`/`npRmpi` terminology onto `crs` or `spcovar`.
- Do not broaden into examples, vignettes, or estimator semantics unless a
  later explicit plan authorizes that scope.

## Acceptance Criteria

The tranche is successful when:

1. every touched package has a pre-edit argument classification ledger;
2. every touched `scale.factor.*` argument entry lists its default clearly;
3. `scale.factor.search.lower = NULL` is explained as inheritance-or-default,
   not as absence of a floor;
4. the four continuous scale-factor controls are contiguous wherever all four
   are public on the page;
5. large pages split overbroad search/kernels/support sections into the
   user-model subsections above;
6. arguments are alphabetical within user-model subsections, and middle
   subsection labels are alphabetical between the leading entry-point section
   and final `Additional Arguments` section;
7. `np` and `npRmpi` wording is semantically equivalent;
8. `crs` and `spcovar` retain or improve existing subsection-header structure
   without adopting irrelevant `np` terminology;
9. every touched page preserves `\usage{}` and argument-name sets unless an
   explicit page note documents a deliberate correction;
10. every touched page has no duplicate argument entries;
11. rendered before/after help output is retained for touched pages;
12. `tools::checkRd()` is clean for every touched `.Rd` file;
13. package build/check documentation gates show no new warnings, notes, or URL
   concerns from this tranche in all touched packages.
