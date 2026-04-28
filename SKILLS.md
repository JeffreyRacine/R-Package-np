# SKILLS.md

## Skill: Rd Mental-Model Documentation Cleanup

Use this skill for `.Rd` documentation updates in `npRmpi`.

1. Work in a tmp campaign worktree, not the live release worktree.
2. Treat documentation cleanup during feature freeze as Rd-only unless explicit
   sign-off authorizes API/formal/MPI runtime changes.
3. For each function family, write the user decision model first:
   - interface/data choice;
   - bandwidth supplied vs selected;
   - bandwidth method/type/search controls;
   - MPI execution mode and lifecycle controls when relevant;
   - kernel and support controls;
   - local-polynomial/degree controls;
   - numerical controls;
   - returned object and follow-up methods.
4. For long `\arguments{}` sections, use `\subsection{Group Name}{}` as an
   empty sibling before top-level `\item{arg}{...}` entries. Real arguments must
   remain top-level `\item{arg}{...}` entries.
5. Do not use `\item{Group Name}{...}` for headings; R will treat it as a
   documented argument and may warn during checks. Do not wrap `\item{}` entries
   inside the second argument of `\subsection{...}{...}`; touched-page
   validation showed that nested form makes `\item{}` parse incorrectly.
6. Keep `\usage{}` synchronized with actual formals. Do not reorder `\usage{}`
   merely for visual tidiness.
7. Do not reorder R function formals in this campaign unless separately
   approved after an explicit positional-use risk audit.
8. On estimator pages, make `bws`/`bandwidth` and `...` explain the implicit
   call to the corresponding `*bw` selector when bandwidths are omitted.
9. Keep estimator pages concise; point to the `*bw` page for the full
   bandwidth/search/kernel/support control surface.
10. Preserve `npRmpi`-specific attach/profile/autodispatch wording when
    mirroring serial documentation structure.
11. Add Gallery pointers only using canonical Gallery wording/URLs already
    present in package materials.
12. Validate each tranche with:
    - `R CMD Rd2txt` for touched pages;
    - `tools::checkRd()` for touched pages;
    - package check from a built tarball/private library;
    - installed smoke calls for estimator pass-through documentation;
    - relevant `nslaves=1`, `nslaves>1`, attach, or profile gates for touched
      MPI routes.
