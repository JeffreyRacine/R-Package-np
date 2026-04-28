# np / npRmpi Release Readiness Plan: Argument/Rd Cleanup First

Date: 2026-04-28

## Goal

Prepare `np` and `npRmpi` for the `0.70-1` release without rushing and without
introducing collateral regressions. The first release-readiness priority is a
careful documentation campaign: function-argument and `.Rd` organization,
argument/function description fidelity, estimator-to-bandwidth discoverability,
and external navigation to the Gallery. The substantive CHANGELOG expansion
should happen after that cleanup so the release notes describe the final public
surface rather than an intermediate one.

Initial CHANGELOG normalization was performed on 2026-04-28 in both repos:
the current top section is `Changes from Version 0.70-0 to 0.70-1
[28-Apr-2026]`. The detailed final release-note expansion remains deferred
until the cleanup and release gates are complete.

## Proof Of Concept: Rd Argument Subsections

An end-to-end proof of concept was completed on 2026-04-28 in a scratch copy of
`np-master`, not in the live package tree:

- scratch package:
  `/Users/jracine/Development/tmp/rd_subsection_poc_np_20260428/np-master-poc`
- modified page:
  `man/np.condensity.bw.Rd`
- change:
  inserted real `\subsection{...}{...}` headings inside `\arguments{}` without
  changing R code, function formals, `.Rd` `\usage{}`, defaults, examples, or
  argument item text.
- rendered text artifact:
  `/Users/jracine/Development/tmp/rd_subsection_poc_np_20260428/np_condensity_bw_poc.txt`
- `tools::checkRd()` artifact:
  `/Users/jracine/Development/tmp/rd_subsection_poc_np_20260428/checkRd_np_condensity_bw_poc.log`
- CRAN-style check artifact:
  `/Users/jracine/Development/tmp/rd_subsection_poc_np_20260428/check_as_cran.log`
- unmodified baseline comparison:
  `/Users/jracine/Development/tmp/rd_subsection_poc_np_20260428_baseline/check_as_cran.log`

Result:

- `R CMD Rd2txt` rendered subsection headings cleanly under `Arguments`.
- `tools::checkRd()` was clean.
- `R CMD check --as-cran` completed with all Rd, usage, contents, manual PDF,
  and manual HTML phases `OK`.
- The proof copy and unmodified baseline both ended with the same package-level
  status, `3 WARNINGs, 1 NOTE`, from existing non-Rd release issues:
  missing `checkbashisms`, vignette/`inst/doc` state, and incoming/date notes.
- No warning, note, or error was introduced by `\subsection{}` inside
  `\arguments{}`.

Conclusion: real `\subsection{}` headings inside `\arguments{}` are an accepted
and useful mechanism for organizing long argument lists. They should be the
preferred first-line documentation enhancement for large public help pages.

## Initial Documentation Surface Observations

A quick structural scan on 2026-04-28 found:

- `np-master` has 57 `.Rd` files.
- `np-npRmpi` has 72 `.Rd` files.
- The longest public pages are argument-heavy and likely to benefit from
  subsection organization:
  - `np.condensity.bw.Rd` with 74 documented items;
  - `np.condistribution.bw.Rd` with 73;
  - `np.regression.bw.Rd` with 61;
  - `np.smoothcoef.bw.Rd` with 54;
  - `np.distribution.bw.Rd` with 50;
  - `np.density.bw.Rd` with 46;
  - `np.singleindex.bw.Rd` with 40.
- Several pages intentionally lack sections such as `\arguments{}` or
  `\details{}` because they are package, data, or helper pages; these should be
  audited semantically rather than "fixed" mechanically.
- No current `.Rd` file matched obvious Gallery terms. If the package now points
  users to a Gallery website, the `.Rd` navigation layer should be reviewed so
  users can discover that living example corpus from appropriate pages.

## Release Principles

1. Preserve estimator behavior exactly unless a failing gate proves a bug that
   must be fixed.
2. Preserve positional-call compatibility for canonical leading arguments.
3. Do not reorder arguments mechanically across the whole package.
4. Do not change defaults, argument names, aliases, dispatch behavior, or
   examples as part of ordering cleanup unless a local mismatch is found and
   validated.
5. Treat `np-master` as the source of truth, then port surgically to
   `np-npRmpi`.
6. Use installed-build validation with private libraries for all positive
   release-readiness claims.
7. Prefer documentation/order clarity over formal-signature churn whenever the
   user-facing gain is the same.
8. Treat every formal-argument reorder as an API change until proven otherwise.
9. Put the reader first: every documentation edit should make a real user
   decision easier, not merely make the source look tidier.
10. Prefer fewer, better headings over dense taxonomy. A help page that becomes
    over-sectioned is not an improvement.
11. Preserve searchability: common argument names, family names, and method
    names should remain easy to find in the rendered text output.

## Critical Review Of The Initial Plan

The initial plan had the right north star but still carried too much implicit
risk for a release candidate.

1. It treated function-formal reordering and `.Rd` argument reordering as a
   paired cleanup. In R, those are not equally risky. `.Rd` ordering is mostly
   presentation; formal ordering can affect positional matching, `match.call()`,
   partial matching edge cases, S3 method compatibility, examples, downstream
   code, and user scripts that rely on positional calls beyond the first few
   arguments.
2. It did not define a strong enough "do not move" rule. Preserving only the
   leading canonical prefix is necessary but not sufficient; some later
   arguments can still be de facto positional in long-lived packages.
3. It did not separate three distinct operations:
   - `.Rd` `\arguments{}` presentation order,
   - `.Rd` `\usage{}` signature order,
   - R function formal order.
   These should not be changed together by default.
4. It allowed cluster-level edits before a baseline artifact existed. For this
   kind of release cleanup, every tranche needs a pre-edit snapshot of formal
   order, usage order, examples, tests, and check status.
5. It did not define rollback criteria sharply enough. A documentation cleanup
   should have a very low tolerance for failed checks, new warnings, changed
   behavior, or changed matching semantics.
6. It underweighted multi-signature `.Rd` pages. Shared help pages can become
   less clear if arguments are alphabetized globally instead of grouped by the
   signatures that actually use them.
7. It did not explicitly reserve mechanical formatting. Reflowing large `.Rd`
   blocks or R signatures while changing order would make review harder and
   increase regression risk.
8. It included performance gates, which are useful for release readiness, but
   argument/Rd cleanup should not be accepted or rejected on performance unless
   it accidentally touches executable logic. Performance gates belong after the
   cleanup and before release, not inside every docs-order tranche.

## Reformulated Strategy

The safer release-grade strategy is:

1. Make `.Rd` argument presentation the primary cleanup surface, using real
   `\subsection{}` headings inside `\arguments{}` for long public pages.
2. Change `.Rd` `\usage{}` only when it must match existing function formals or
   when the function formal order has deliberately passed the formal-risk gate.
3. Change R function formal order only when all of these are true:
   - the function is public-facing enough that the clearer order matters;
   - the arguments to be moved are demonstrably option-like rather than
     canonical positional arguments;
   - no S3 generic signature, internal dispatcher, `match.call()` consumer, or
     test fixture depends on the old order;
   - focused positional-call probes pass before and after;
   - the change can be reviewed in a small diff.
4. Leave internal helpers and low-level C/R bridge wrappers alone unless their
   `.Rd` page is user-facing and can be clarified without changing formals.
5. Use one family/page cluster per tranche, with a hard rollback if the tranche
   creates any check warning/error or any behavioral/positional mismatch.
6. Defer the final substantive CHANGELOG expansion until the cleanup outcome is
   known and validated.

## Final Risk Refinement For Feature Freeze

Because the packages are at feature freeze and the goal is release hardening,
the default campaign posture is now stricter than ordinary documentation
cleanup:

1. No R formal-argument reordering in the first implementation pass.
2. No `.Rd` `\usage{}` reordering in the first implementation pass unless it is
   needed to repair an existing documented/code mismatch.
3. No default changes, no alias changes, no signature changes, no S3 method
   signature changes, and no changes to runtime dispatch.
4. `.Rd` `\arguments{}` subsection organization, argument text, `\details{}`,
   `\seealso{}`, and lightweight examples are the allowed first-pass surfaces.
5. Any proposed R formal-order change must be spun out into a separate
   post-inventory decision with explicit sign-off, because even cosmetic formal
   movement can break downstream positional code.
6. The first live documentation edit should be a pilot family, not a broad
   sweep. Recommended pilot: `npcdens`/`npcdensbw`, because it exercises the
   hardest documentation problems: estimator-to-`*bw` pass-through, bounded
   support, local-polynomial controls, scale-factor controls, quadrature, and
   NOMAD/Powell search.
7. If the pilot produces any new check warning/error, any confusing rendered
   help output, or any smoke-test mismatch, revert the pilot tranche before
   touching another family.
8. Commit or checkpoint only after a tranche passes its local gates. Do not
   accumulate multiple unrelated family edits before validation.
9. If a tranche fails, revert that tranche in the tmp worktree before starting a
   different family; do not layer cleanup on top of a known-bad docs state.

Branching state for this campaign:

- `np-master`: use a dedicated documentation branch such as
  `codex/rd-doc-mental-model-np`.
- `np-npRmpi`: use the paired documentation branch
  `codex/rd-doc-mental-model`.
- Keep the branches separate because the repos share a worktree-style topology
  and may not be able to use the same branch name simultaneously.

## Rd Subsection Style Contract

Use `\subsection{Group Name}{...}` inside `\arguments{}` for long public pages
whose argument lists are hard to scan. Do not fake headings with
`\item{Group Name}{...}` because that creates documented non-arguments and
triggers `R CMD check` usage warnings.

Subsection text should be short and user-facing. It should say what kind of
decisions the group controls, not restate every argument.

Recommended group names should be stable across related pages when applicable:

1. `Data And Bandwidth Inputs`
2. `Formula Interface`
3. `Bandwidth Selection Controls`
4. `Continuous Scale-Factor Search Initialization`
5. `Discrete And Categorical Search Initialization`
6. `Kernel Type Controls`
7. `Continuous Kernel Support Controls`
8. `Local-Polynomial And Degree Search Controls`
9. `Least-Squares Quadrature Controls`
10. `Numerical Search And Tolerance Controls`
11. `Plot And Interval Controls`
12. `Returned Values`
13. `Compatibility And Additional Arguments`

Use only the groups that make sense for the page. Prefer five or fewer groups
for moderate pages; use more only for very large pages such as bandwidth
selectors with kernel, search, support, and local-polynomial controls.

Reader-experience rules:

1. A subsection heading should answer "what decision am I making here?"
2. A subsection should usually contain at least two real arguments; avoid
   one-argument headings unless the argument is genuinely central and complex.
3. Keep high-use arguments near where users expect them, even if a different
   taxonomy would be more internally elegant.
4. Do not bury `formula`, `data`, `bws`, `xdat`, `ydat`, or `...`.
5. Rendered `Rd2txt` output is decisive for readability; source organization is
   secondary.

## Documentation Style Ledger

Before broad editing, create and maintain a campaign-local style ledger:

`issue_notes/rd_documentation_style_ledger_20260428.md`

The ledger should be short, reusable, and updated only when a wording or
structure pattern is actually approved by a passing tranche. It should record:

1. approved subsection names and their intended meaning;
2. approved estimator-to-`*bw` pass-through wording;
3. approved `...` wording patterns by page type;
4. approved `bws`/`bandwidth` wording patterns;
5. approved Gallery wording and URL once confirmed from existing vignettes;
6. approved `\seealso{}` ordering patterns;
7. examples of wording that should not be used, especially stale control names
   or fake `\item{}` headings.

Do not turn the ledger into a status tracker. It is a style contract for future
documentation edits.

## Estimator-to-Bandwidth Documentation Contract

Many high-level estimator functions (`npreg`, `npudens`, `npcdens`,
`npudist`, `npcdist`, `npplreg`, `npindex`, `npscoef`, and related formula
front ends) accept a compact visible argument set, but when called without an
already-computed bandwidth object they invoke their `*bw` counterpart and pass
bandwidth-selection controls through `...`. This is powerful but currently too
easy for users to miss.

During each estimator/Rd cleanup tranche, add a small, consistent bridge from
the estimator page to its bandwidth-selector page.

Required estimator-page pattern:

1. In the `bws` or `bandwidth` argument item, explicitly state that a supplied
   bandwidth object is used directly, and that if it is omitted the estimator
   calls the corresponding `*bw` function.
2. In the `...` argument item, explicitly state that when bandwidths are
   omitted, bandwidth-selection controls supplied through `...` are passed to
   the corresponding `*bw` function.
3. Add a short `\subsection{Bandwidth Selection Arguments}{...}` under
   `\arguments{}` or in `\details{}` when the pass-through behavior is otherwise
   easy to overlook.
4. Point users to the `*bw` help page for the full bandwidth, kernel, support,
   search, local-polynomial, quadrature, and scale-factor control surface.
5. Ensure `\seealso{}` prominently links to the corresponding `*bw` page.
6. Add or update one compact example where useful, showing a high-level
   estimator receiving a bandwidth-control argument directly, for example:
   `npcdens(y ~ x, data = dat, bwmethod = "cv.ls",
   scale.factor.search.lower = 1e-6)`.

Do not copy the full `*bw` argument universe into each estimator page. The
estimator page should explain the pass-through contract and identify the `*bw`
page as the authoritative reference. This keeps estimator docs readable while
making the hidden bandwidth-control surface discoverable.

Validation for every estimator page updated under this contract:

- render/check the estimator page and its `*bw` counterpart;
- run at least one installed-build smoke call where a bandwidth-selection
  argument is supplied to the estimator without precomputed bandwidths;
- confirm the corresponding explicit `*bw` call still gives matching bandwidth
  behavior for the chosen smoke case.

## Documentation Fidelity And Navigation Contract

This campaign is not only an argument-order campaign. Each touched `.Rd` page
should be checked against the current code so the documentation faithfully
describes what the function does, what each argument controls, and where users
should go for examples or deeper guidance.

For every touched page, audit:

1. Argument existence:
   - every documented argument appears in at least one `\usage{}` signature, or
     has an intentional explanation for multi-signature/help-page structure;
   - every user-facing formal argument is documented, including arguments
     accepted through public methods and meaningful `...` pass-throughs.
2. Default fidelity:
   - documented defaults match the R formals or the documented runtime default
     resolution;
   - `NULL` defaults are explained when they resolve to computed defaults.
3. Behavioral fidelity:
   - descriptions match current code paths, not historical behavior;
   - removed/deprecated controls are described only where they are still
     intentionally trapped with errors or compatibility guidance.
4. Pass-through fidelity:
   - `...` is described concretely, especially estimator-to-`*bw` pass-through
     and plot/helper pass-through behavior.
5. Examples:
   - examples are small enough for checks, still run, and demonstrate current
     recommended idioms;
   - heavier or richer examples should point to the Gallery rather than
     overloading `.Rd` examples.
6. Cross-links:
   - `\seealso{}` links to the corresponding estimator, `*bw`, plot/predict,
     kernel, and options pages where useful;
   - package-level or family-level pages link to the Gallery when appropriate.
7. Terminology:
   - use consistent names for scale-factor controls, bandwidth selectors,
     bounded-kernel support, local-polynomial controls, NOMAD/Powell search,
     and MPI modes;
   - avoid stale names such as superseded `cfac.init` terminology except where
     documenting migration from old names.
8. np/npRmpi parity:
   - shared estimator documentation should be semantically equivalent;
   - `npRmpi` pages must preserve MPI-specific lifecycle, attach/profile, and
     autodispatch guidance rather than copying serial text blindly.

Reader QA rubric for each touched public page:

1. Can a user identify the main call pattern within 10 seconds?
2. Can a user tell whether they should call the estimator or the `*bw`
   counterpart?
3. Can a user find where bandwidth/search controls are documented?
4. Can a user distinguish estimator controls from bandwidth-selection controls?
5. Can a user understand what `...` does on this page?
6. Can a user find bounded-kernel, local-polynomial, NOMAD/Powell, and
   scale-factor controls without reading the whole page?
7. Can a user find the next place to go: estimator, bandwidth selector,
   kernel/options page, vignette, or Gallery?

If the answer to any applicable question is "no", the page is not done even if
`R CMD check` passes.

Gallery integration:

- Do not turn `.Rd` pages into long tutorials when the Gallery is the better
  place for evolving examples.
- Extract the canonical Gallery URL/text from the existing vignettes; do not
  invent a new wording or destination.
- Add concise Gallery pointers from package-level, getting-started-adjacent,
  and high-use family pages only after the canonical Gallery wording is
  confirmed from the current package materials.
- Prefer wording such as "For extended examples and visual workflows, see ..."
  in `\details{}` or `\seealso{}`, keeping CRAN examples lightweight.

## Argument Ordering Contract

Function formals and `.Rd` argument sections should use a consistent grouping
scheme while preserving the historically important leading call shape.

### Leading Positional Prefix

Keep these first when they are already first or clearly canonical:

- formula/data interfaces: `formula`, `data`, `subset`, `na.action`
- bandwidth-object interfaces: `bws`, `bandwidth`, `object`, `x`
- low-level data interfaces: `xdat`, `ydat`, `txdat`, `tydat`, `exdat`, `eydat`
- prediction/evaluation interfaces: fitted object first, then `newdata` or
  route-specific evaluation data
- S3 methods: required generic signature order such as `x, ...` or `object, ...`

The prefix is a compatibility boundary. Reordering may begin only after the
canonical prefix for that function/page is fixed.

### Formal-Order Risk Gate

Before changing a function's formal order, record and approve:

- whether the function is exported;
- whether it is an S3 generic or S3 method;
- whether its examples or tests use positional arguments beyond the canonical
  prefix;
- whether implementation code uses `match.call()`, `missing()`,
  positional forwarding through `...`, or calls to sibling functions that rely
  on the current order;
- whether aliases or `.Rd` pages share the same signature;
- whether downstream scripts under `/Users/jracine/Development/R`,
  `demo/`, `tests/`, and `benchmarks/` call the function positionally.

If any item is uncertain, do not change formal order in that tranche. Reorder
only `.Rd` `\arguments{}` presentation and document the reason in the audit
table.

### Mental-Model Groups After The Prefix

Within each group, sort arguments alphabetically only when doing so improves
scannability. User mental model wins over strict alphabetization. For example,
bound controls should usually read as `*kerbound`, `*kerlb`, `*kerub`, and
degree-search controls should usually keep natural min/max/start/restart/verify
clusters together even if a strict alphabetical order would separate them.

1. Bandwidth/search controls:
   - `bwmethod`, `bwscaling`, `bwtype`, `nmulti`, `remin`, `itmax`,
     `search.engine`, `nomad`, `scale.factor.*`
2. Kernel controls:
   - continuous, unordered, and ordered kernel type arguments and their
     X/Y-specific variants
3. Bound/support controls:
   - `*kerbound`, `*kerlb`, `*kerub`, and related finite-support arguments
4. Local-polynomial controls:
   - `regtype`, `basis`, `degree`, `degree.*`, `bernstein.*`
5. Evaluation/grid controls:
   - evaluation data, grids, quadrature controls, proper-density controls
6. Plot/bootstrap/interval controls:
   - `plot.*`, `gradients`, `errors`, `boot.*`, `B`, confidence controls
7. Computation/progress controls:
   - `memfac`, `ftol`, `tol`, `small`, `large`, progress/message options
8. Miscellaneous compatibility controls:
   - legacy aliases retained for documented use, then `...` last when present

## Tranche 0: Branch And Baseline Safety

1. Work only on the dedicated documentation branches.
2. Record dirty-tree status before each tranche.
3. Build private-library baselines before the first live `.Rd` edit:
   - `R CMD check --no-manual --ignore-vignettes` for fast local baseline;
   - retain current known `--as-cran` status separately so pre-existing release
     warnings are not confused with documentation regressions.
4. Save all logs and generated inventories under
   `/Users/jracine/Development/tmp/release_rd_docs_YYYYMMDD`.
5. Confirm no simulation/manuscript/campaign files are touched.
6. Create the initial documentation style ledger with provisional patterns
   marked as provisional until the pilot tranche validates them.

Acceptance:

- Both repos are on dedicated documentation branches.
- Baseline status is recorded.
- No package source edits have been made beyond existing CHANGELOG/plan work.
- The style ledger exists and clearly distinguishes provisional from validated
  patterns.

## Tranche 1: Inventory Ledger, Baseline, And Risk Map

1. Generate a formal inventory for both repos:
   - exported functions and S3 methods;
   - function formal order from `R/`;
   - `.Rd` `\usage{}` signatures;
   - `.Rd` `\arguments{}` item order.
2. Classify pages into risk tiers:
   - Tier A: high-use public bandwidth/estimator/plot pages.
   - Tier B: helper/public utility pages.
   - Tier C: internal or rarely used documented pages.
3. Identify pages with multiple aliases/usages where one shared `\arguments{}`
   block must serve several signatures.
4. Search examples, tests, demos, benchmarks, and `/Users/jracine/Development/R`
   for positional use of Tier A functions.
5. Search R code for `match.call()`, `missing()`, and positional forwarding in
   Tier A families.
6. Produce a no-edit audit table showing:
   - current formal order;
   - current `.Rd` usage order;
   - current `.Rd` argument order;
   - alias list and family membership;
   - argument count and subsection need;
   - whether `...` exists and how it is currently described;
   - whether a `bws`/`bandwidth` object exists and how it is described;
   - likely `*bw` counterpart, if any;
   - current `\seealso{}` targets;
   - current Gallery/vignette/example navigation;
   - proposed `.Rd` argument grouping;
   - proposed user decision model summary;
   - documentation fidelity issues;
   - missing or weak `...`/pass-through descriptions;
   - missing or weak `\seealso{}`/Gallery navigation;
   - formal-order risk rating;
   - explicit "formals may move" or "Rd-only" recommendation.
7. Run a baseline package check for each repo from a private library before the
   first edit tranche.
8. Generate a "no-touch" list:
   - low-level/internal helper pages where cleanup value is lower than review
     risk;
   - pages whose current structure is intentionally sparse, such as data pages;
   - any page requiring statistical wording decisions rather than mechanical
     documentation cleanup.

Acceptance:

- No code edits.
- Audit table reviewed against examples and tests.
- The inventory ledger is saved as CSV/TSV plus a short Markdown summary under
  `/Users/jracine/Development/tmp/release_rd_docs_YYYYMMDD`.
- Inventory scripts and baseline logs retained under `/Users/jracine/Development/tmp`.
- Baseline package checks either pass or have documented pre-existing issues
  unrelated to the cleanup.
- Each Tier A family has a draft user decision model before any `.Rd` edit.

## Tranche 1A: Family Mental-Model Prewrite

Before editing a family, write a compact mental-model note for that family in
the campaign worktree, for example:

`issue_notes/rd_mental_model_npcdens_20260428.md`

Each note should answer:

1. What is the user's sequence of decisions for this function family?
2. Which arguments belong to each decision?
3. Which page is authoritative for bandwidth/search controls?
4. Which estimator pages need pass-through wording?
5. Which arguments have computed defaults or runtime-resolved `NULL` behavior?
6. Which examples should remain in `.Rd`, and which belong in the Gallery?
7. Which `npRmpi` additions are MPI-specific rather than shared with `np`?

Acceptance:

- No `.Rd` edits before the family mental-model note exists.
- The note is concise enough to guide edits without becoming a tutorial.
- The note identifies any statistical wording decisions that need Jeffrey's
  review before implementation.

## Tranche 2: Pilot Family In np, Rd-Only

Start with one `np-master` pilot family, preferably `npcdens`/`npcdensbw`.
This is deliberately Rd-only: no R formal changes and no `.Rd` `\usage{}`
reordering unless repairing a pre-existing mismatch.

Pilot tasks:

1. Create the family mental-model note.
2. Add mental-model `\subsection{}` organization to the `*bw` page.
3. Add estimator-to-`*bw` pass-through wording to the estimator page.
4. Reconcile touched argument/default text with current code.
5. Add or improve `\seealso{}` links.
6. Add Gallery navigation only using canonical text/URL from existing vignettes.
7. Render before/after `R CMD Rd2txt` artifacts for touched pages.
8. Run `tools::checkRd()` on touched pages.
9. Build/check from a tarball with private library.
10. Run installed smoke calls covering:
   - ordinary estimator call;
   - estimator call with a bandwidth-selection argument passed through `...`;
   - explicit `*bw` call for the same bandwidth-control smoke case.
11. Update the style ledger only for wording/structure patterns proven by the
    pilot.
12. Perform a rendered-help reader review using the Reader QA rubric.
13. Record the before/after readability conclusion in the pilot artifact
    summary.

Pilot acceptance:

- No R files changed.
- No `\usage{}` changes unless justified by a pre-existing mismatch.
- `R CMD Rd2txt` output is clearer and not cluttered.
- Before/after rendered help artifacts are retained.
- `tools::checkRd()` is clean for touched pages.
- `R CMD check --no-manual --ignore-vignettes` has no new warning/error/note.
- Installed smoke calls pass and match the intended `*bw` pass-through
  behavior.
- The rendered help tells a coherent user decision story.
- The style ledger marks pilot-proven wording as validated.
- The pilot artifact summary states why the new help is easier for a user to
  navigate, not only that checks passed.

## Tranche 3: Remaining High-Use np Pages, Rd-First With Subsections

After the pilot passes, continue in `np-master` one family or page cluster at a
time:

1. `npreg*`
2. `npudens*`
3. `npcdens*`
4. `npudist*`
5. `npcdist*`
6. `npscoef*`
7. `npindex*`
8. `npplreg*`
9. plot/predict/summary pages adjacent to touched families

For each page/function:

1. Record pre-edit formal order and `.Rd` order.
2. Decide the positional prefix.
3. Reorder `.Rd` `\arguments{}` by user mental-model groups, using
   `\subsection{}` headings when the page is long enough to benefit.
4. For estimator pages, add the estimator-to-bandwidth documentation bridge:
   clear `bws`/`...` pass-through wording, a `*bw` cross-reference, and one
   compact example where useful.
5. Audit and repair documentation fidelity for touched pages:
   argument descriptions, defaults, `...` pass-throughs, examples, links, and
   current-code behavior.
6. Add Gallery pointers where appropriate once canonical wording/URL is agreed.
7. Leave `.Rd` `\usage{}` matching current function formals unless repairing a
   pre-existing mismatch.
8. Do not reorder R formals in this campaign unless a separate sign-off
   explicitly authorizes a post-inventory formal-order tranche.
9. Avoid unrelated line wrapping or formatting churn.
10. Run focused examples/tests for that family.
11. Run the Reader QA rubric on rendered help for the edited public pages.
12. Run `R CMD check --no-manual --ignore-vignettes` after each cluster.
13. Commit/checkpoint the cluster only after the rendered help review and
    package checks pass.

Acceptance:

- No changed defaults.
- No changed argument names.
- No changed S3 generic/method required signatures.
- No new `R CMD check` warnings/errors.
- Focused positional smoke tests pass for canonical calls.
- `R CMD Rd2txt` confirms subsection headings render cleanly for at least one
  page in each edited family.
- Help pages render and usage blocks remain consistent.
- Estimator pages that compute bandwidths clearly identify the corresponding
  `*bw` page as the complete bandwidth-control reference.
- At least one estimator-with-pass-through smoke call passes per edited family.
- Touched argument descriptions and defaults are checked against current code.
- Touched pages have intentional `\seealso{}`/Gallery navigation.
- Diff is small enough to review function by function.
- Rendered help passes the Reader QA rubric for the edited public pages.
- The tranche is independently revertible.

## Tranche 4: Mirror Into npRmpi Surgically

After an `np-master` cluster passes:

1. Port only the corresponding argument/Rd ordering changes.
2. Preserve `npRmpi`-specific MPI arguments, examples, and lifecycle wording.
3. Do not copy whole R or `.Rd` files between repos.
4. Preserve `npRmpi`-specific formal order when MPI lifecycle or autodispatch
   arguments make the serial order inappropriate.
5. Validate each mirrored cluster in:
   - serial/no-slave mode where applicable;
   - attach/autodispatch mode for touched MPI routes;
   - profile/manual-broadcast mode for at least one representative route per
     touched family.
6. Update the style ledger only when `npRmpi` requires a distinct MPI wording
   pattern.

Acceptance:

- `npRmpi` remains standalone, with no runtime `np::` bridge calls in package
  internals.
- `R CMD check --no-manual --ignore-vignettes` passes.
- No MPI worker or `mpiexec` leftovers after validation.

## Tranche 5: Remaining Public And Helper Pages

Handle Tier B/C pages after the core family surface is stable.

Rules:

- Prefer `.Rd` ordering-only changes for helper pages unless function formal
  order is clearly public-facing and safe to adjust.
- Leave low-level/internal signatures alone if reordering would create more
  risk than clarity.
- If a page documents multiple signatures with partially disjoint arguments,
  group by shared leading prefix first, then by option families.
- For internal helpers documented only for maintainers, accept "leave as-is" as
  a valid outcome when the cleanup value is smaller than the review risk.
- For package/data/helper pages without `\arguments{}` or `\details{}`, decide
  case by case whether adding a section improves user understanding; do not add
  empty boilerplate merely for uniformity.
- Bring `\seealso{}` coverage up to an intentional standard for public pages,
  especially where users need to discover bandwidth selectors, estimator
  counterparts, kernel references, options, or Gallery examples.

Acceptance:

- Full package checks pass in both repos.
- No usage/documentation mismatch.
- No example drift.
- Documentation omissions are either fixed or explicitly classified as
  intentional/no-action in the audit table.

## Tranche 6: Substantive CHANGELOG Expansion

After argument/Rd cleanup and validation:

1. Update the top `0.70-1` sections in both CHANGELOG files.
2. Add major changes since the prior update above the existing top bullets.
3. Keep bullets user-facing and grouped by release theme rather than commit log.
4. Suggested high-level groups:
   - scale-factor search controls and lower-bound policy;
   - bounded conditional-density CV-LS quadrature and NOMAD/Powell repairs;
   - `scale.factor.*` naming modernization;
   - local-linear raw-basis failure guidance;
   - bounded convolution extinct-code excision;
   - large-bandwidth fast gates and summary accounting;
   - MPI lifecycle, streaming, attach/profile validation, and CRAN unload race
     remediation for `npRmpi`;
   - documentation/API cleanup from the argument/Rd campaign.
5. Recheck spelling, dates, and repo-specific wording.

Acceptance:

- CHANGELOGs mention only validated changes.
- `np` and `npRmpi` entries are consistent but not blindly identical.
- Final top header remains `Changes from Version 0.70-0 to 0.70-1
  [28-Apr-2026]` unless the release date changes.

## How To Put The Packages Through Their Paces

### Gate 1: Structural/API Audit

- Formal/Rd inventory diff before and after cleanup.
- Documentation inventory ledger before/after each family tranche.
- Family mental-model note exists before edits.
- `tools::checkRd()` or package check Rd phase via `R CMD check`.
- `R CMD Rd2txt` spot checks for representative subsectioned pages, especially
  the largest bandwidth selector pages.
- Before/after rendered help artifacts for the pilot and any unusually large
  page.
- Documentation bridge audit: each estimator page with implicit bandwidth
  computation must clearly link its `...` bandwidth-control pass-through to the
  corresponding `*bw` page.
- Documentation fidelity audit: touched pages must have argument/default
  descriptions reconciled with current code, intentional examples, and useful
  cross-links/Gallery navigation.
- Reader QA rubric must pass on rendered help for touched public pages.
- Positional smoke tests for canonical leading calls:
  `npreg(y ~ x, data=...)`, `npreg(bws=...)`, `npudens(~ x, data=...)`,
  `npcdens(y ~ x, data=...)`, `npudist(~ x, data=...)`,
  `npcdist(y ~ x, data=...)`, `npplreg`, `npindex`, `npscoef`.
- Extra positional smoke tests for any function whose formal order was changed,
  including at least one call that uses positional arguments through the last
  moved argument.

### Gate 2: Installed-Build Functional Sentinels

For `np`:

- focused testthat families for each touched route;
- simple sentinel scripts from `demo/`, `tests/`, and `/Users/jracine/Development/R`;
- representative heavier `npcdens` bounded CV-LS/NOMAD and semiparametric
  search-control sentinels.

For `npRmpi`:

- the same family sentinels where mirrored;
- `nslaves=1` and `nslaves=2` runs for representative `npreg`, `npcdens`,
  `npcdist`, `npindex`, and `npscoef` paths;
- attach/autodispatch and profile/manual-broadcast gates for touched MPI routes.

### Gate 3: Performance/Scaling Smoke

- Run after argument/Rd cleanup gates pass, not as a substitute for them.
- Do not broad-brute-force first.
- Use structural comparability to select suspicious routes.
- For suspicious routes, run paired pre/post timing with private installed
  libraries and saved raw logs.
- For `npRmpi`, include at least one scaling probe where serial `np`,
  `npRmpi nslaves=1`, and `npRmpi nslaves=2` are compared on a problem large
  enough that MPI overhead should not dominate.

### Gate 4: CRAN And Platform Checks

- `R CMD check --as-cran` tarball-first for `np`.
- `R CMD check --as-cran` tarball-first for `npRmpi` where local MPI allows;
  otherwise retain local `--no-manual --ignore-vignettes` plus targeted MPI
  lifecycle gates.
- win-builder / rhub-equivalent checks as available.
- Explicit load/unload subprocess checks for `npRmpi`.

### Gate 5: Final Release Rehearsal

- Build final tarballs.
- Install into clean private libraries.
- Run named sentinel bundle from installed packages only.
- Verify `sessionInfo()`, package versions, CHANGELOG dates, and no dirty
  worktrees.
- Keep raw logs under `/Users/jracine/Development/tmp/release_readiness_YYYYMMDD`.
- Perform one final rendered-help spot review across the largest bandwidth
  selector pages and the estimator pages that point to them.
- Confirm the final documentation campaign changed no R runtime files unless a
  separately approved tranche explicitly did so.

## Immediate Next Step

Create the campaign style ledger, then create the formal/Rd inventory script and
audit table for `np-master`, including the formal-order risk gate fields and
documentation-fidelity columns. Do not make any ordering edits until the Tier A
inventory and the first family mental-model note exist. Do not touch `npRmpi`
until the corresponding `np-master` tranche passes.
