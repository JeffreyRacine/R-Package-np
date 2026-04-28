# npcdens Rd Mental Model Pilot For npRmpi

Date: 2026-04-28

Scope: documentation-only pilot for `npcdens()` and `npcdensbw()` in `npRmpi`,
ported surgically from the validated `np` pilot while preserving MPI-specific
startup/lifecycle wording.

## User Decision Sequence

1. Initialize and manage MPI execution as appropriate for the session
   (`npRmpi.init()`, attach/profile demos, and `npRmpi.quit()`).
2. Choose the interface: formula, data-frame, or precomputed bandwidth.
3. Decide whether bandwidths are supplied or computed. If `bws` is omitted,
   `npcdens()` calls `npcdensbw()` and passes bandwidth-selection controls
   through `...`.
4. Choose the bandwidth criterion: `cv.ml`, `cv.ls`, or `normal-reference`.
5. Choose kernel/support handling. Continuous `xdat` controls begin with
   `cxker*`, continuous `ydat` controls begin with `cyker*`, and categorical
   controls begin with `u*ker*` or `o*ker*`.
6. Tune search only when needed. `scale.factor.init`,
   `scale.factor.init.lower`, and `scale.factor.init.upper` govern
   continuous search starts; `scale.factor.search.lower` governs the hard
   admissibility floor for continuous fixed-bandwidth search candidates.
7. Use local-polynomial/NOMAD controls only when the polynomial-adaptive route
   is intended.

## Documentation Contract

`npcdensbw()` is the authoritative help page for bandwidth-selection, kernel,
support, local-polynomial, quadrature, NOMAD, and scale-factor search controls.
`npcdens()` should keep its argument list focused on estimation inputs and
state that bandwidth-selection controls are accepted through `...` when
bandwidths are omitted.

The `.Rd` pilot must not reorder R function formals, alter `\usage{}`
signatures, change defaults, or change MPI execution semantics. Argument group
headings must use the Rd-valid empty sibling form
`\subsection{Group Name}{}` followed by top-level `\item{arg}{...}` entries.
