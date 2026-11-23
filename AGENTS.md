# AGENTS.md

Use `/Users/jracine/Development/AGENTS.md` as the canonical workflow and benchmarking policy.
Follow its sweep safety rules and incorporated protocols. For release qualification, candidate replacement and systemwide installation verification, use the workspace task guide at `/Users/jracine/Development/SKILLS.md`.

Repo-specific note:
- Merge (do not overwrite) R-layer and man-page files when porting from `np-master`.
- Any new or ported benchmark/function must follow MPI protocol from `/Users/jracine/Development/AGENTS.md`: start/stop slaves, broadcast commands/options/data explicitly, support passable `nslaves`/`rslaves`, and keep parity with serial seed/DGP/options.
- For release timing alerts, use the workspace's bounded, same-host paired comparison. Report seconds and percent for the affected route, tree setting and worker count. Where MPI shows a jump, separate the changed computation from communication before attributing a slowdown.
- For optimization claims, follow the paired center, tail, scaling and sample-size rules in `/Users/jracine/Development/protocols/PERFORMANCE_RULES.md`. Start MPI screening with one worker and include multiple workers for a scaling claim. Spawn sessions require at least one worker.
- Keep shared regression summary telemetry in sync with `np-master`: `num.feval`, `num.feval.fast`.
- Keep `man/np.kernels.Rd` and kernel/options/plot cross-links aligned, while preserving package-specific `npRmpi` docs text.
- Runtime contract: `npRmpi` must never call serial `np` estimator/plot/bootstrap/helper code paths and must not silently fall back to serial execution when MPI mode is selected.
- Release-hardening reminder: the 2026-05-01 CRAN pretest caught a real
  `np` reverse-dependency regression in bandwidth-dispatch compatibility after
  local checks passed. Matching compatibility repairs were mirrored to
  `npRmpi`. Any change to formula wrappers, `...` forwarding,
  bandwidth-constructor dispatch, call materialization, compatibility helpers,
  or MPI materialization of those surfaces must be validated through the
  revdep-aware release gate in `/Users/jracine/Development/release_protocol`
  plus the relevant `npRmpi` installed/MPI smoke gates.
- Native-code release-hardening reminder: for changes touching `src/`,
  registered native interfaces, MPI/native bridge payloads, or
  `.C`/`.Call`/`.Fortran` payload lifetimes, run the release gate with
  `RUN_RCHK=1` when container/MPI infrastructure is available. `RUN_RCHK=auto`
  is acceptable for ordinary rehearsal only if the resulting summary records
  either PASS or a precise SKIP reason, including missing MPI tooling in the
  rchk container when applicable.
- Demo execution contract (`np-npRmpi/demo`):
  - treat `demo/makefile` as launcher source of truth for attach/profile;
  - run from `demo` subdirectories (`serial`, `n_*_attach`, `n_*_profile`) using `make -f ../makefile ...`;
  - in-repo profile path should resolve to `../inst/Rprofile`; for copied demo folders outside repo, set `RPROFILE` explicitly to an absolute profile path (typically `system.file("Rprofile", package="npRmpi")`).
  - for profile/manual-broadcast launch, provide exactly one profile source (`.Rprofile` in cwd or `R_PROFILE_USER`); never set both `R_PROFILE_USER` and `R_PROFILE` to the same file in a single `mpiexec` launch.
  - `R CMD BATCH --no-save` is valid for profile mode when the single-source profile rule is followed.
  - profile startup guard in `inst/Rprofile` is expected to fail-fast on dual-source profile wiring with a remediation message.
- Demo failure recovery contract:
  - reproduce first with `NP_DEMO_N=100`, `NP=2`, and `DEMOS='<single demo>'`;
  - validate effective profile command with `make -n run-profile` and confirm `R_PROFILE_USER` path;
  - for suspected hangs, monitor worker split with `ps -p '<pid1,pid2>' -o 'pid,stat,etime,%cpu,command'`;
  - inspect progress and process ownership before treating an idle rank as a hang. For a confirmed hang, stop the owning launcher and clean up only its verified processes, following the workspace MPI rules.

## Rd Documentation Campaign Policy (2026-04-28)

During release hardening, `.Rd` documentation cleanup is a user-facing
readability campaign, not an API or MPI-runtime refactor.

1. Organize long `\arguments{}` sections by the user's decision model for the
   function, not by historical/internal argument order.
2. Use real, non-empty sibling subheaders such as
   `\subsection{Group Name}{Short orientation sentence.}` before top-level
   `\item{arg}{...}` entries when long `\arguments{}` sections need group
   headings. Do not wrap `\item{}` entries inside `\subsection{...}{...}`, do
   not leave subheaders empty, and do not fake headings with
   `\item{Group Name}{...}`.
3. First-pass documentation cleanup is Rd-only:
   - no R formal argument reordering,
   - no `.Rd` `\usage{}` reordering unless fixing an existing mismatch,
   - no defaults, aliases, signatures, S3 signatures, MPI lifecycle semantics,
     autodispatch behavior, or runtime dispatch changes.
4. Estimator pages that compute bandwidths implicitly must explain that
   bandwidth-selection controls in `...` pass through to the corresponding
   `*bw` function, and must point users to the `*bw` help page as the complete
   bandwidth-control reference.
5. Preserve `npRmpi`-specific lifecycle, attach/profile, autodispatch, and
   slave-cleanup guidance when mirroring serial `np` documentation structure.
6. For long public pages, add a concise argument-group guide in `\details{}`
   when the plain `\arguments{}` list is too long to reveal the user decision
   model by item order alone.
7. Use the standard heading `Search Initialization, Kernels, And Support` for
   the shared bandwidth-search group; all continuous and categorical kernel
   arguments for that page belong in this group. Keep local-polynomial and
   NOMAD controls in a separate group when those controls exist.
8. For touched pages, verify argument descriptions, defaults, `...`
   pass-throughs, examples, and `\seealso{}` links against current code.
9. Gallery links should use canonical wording/URLs from existing package
   materials; do not invent new Gallery destinations while editing `.Rd` pages.
10. Validate documentation edits with `R CMD Rd2txt`, `tools::checkRd()`, package
   checks from a tarball/private library, installed smoke calls, and relevant
   MPI attach/profile or `nslaves` gates for touched routes.
11. Any proposed R formal-order change is a separate post-inventory decision
   requiring explicit sign-off because it can break downstream positional code.

## Canonical Implementation Directive (2026-03-05)

This repository follows a strict canonical execution rule:

1. Preserve the statistical contract of each method across its qualified computational paths. The workspace engineering context documents the supported dense and tree paths.
2. Unsupported configurations must fail fast with explicit `stop(...)` diagnostics.
3. No silent remap/coercion of user-selected options (for example `bwmethod`, `regtype`, kernels, `cv.iterate`, or bounds transforms).
4. No hidden alternate execution paths for the same method semantics.
5. All fit-defining options (for example `degree`, `basis`, `bernstein.basis`, kernels, and bounds) must be propagated and used by the canonical path.
6. The default is `np.tree = "auto"`. Preserve the explicit `TRUE` and `FALSE` settings, including the documented eligibility rules. For selector changes, verify numerical parity and the computation actually selected in each affected mode.
7. Remove or reject legacy/debug compatibility branches that add redundant runtime overhead once canonical behavior exists.
