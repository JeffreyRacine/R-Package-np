# npcdens Rd Mental Model Pilot

Date: 2026-04-28

Scope: documentation-only pilot for `npcdens()` and `npcdensbw()` in `np`.

## User Decision Sequence

1. Choose the interface: formula, data-frame, or precomputed bandwidth.
2. Decide whether bandwidths are supplied or computed. If `bws` is omitted,
   `npcdens()` calls `npcdensbw()` and passes bandwidth-selection controls
   through `...`.
3. Choose the bandwidth criterion: `cv.ml`, `cv.ls`, or `normal-reference`.
4. Choose kernel/support handling. Continuous `xdat` controls begin with
   `cxker*`, continuous `ydat` controls begin with `cyker*`, and categorical
   controls begin with `u*ker*` or `o*ker*`.
5. Choose the local method: `lc`, legacy canonical `ll`, or explicit `lp`.
6. Tune search only when needed. `scale.factor.init`,
   `scale.factor.init.lower`, and `scale.factor.init.upper` govern
   continuous search starts; `scale.factor.search.lower` governs the hard
   admissibility floor for continuous fixed-bandwidth search candidates.
7. Use bounded `cv.ls` quadrature controls only when bounded continuous
   response support and the `cv.ls` route make them relevant.

## Documentation Contract

`npcdensbw()` is the authoritative help page for bandwidth-selection, kernel,
support, local-polynomial, quadrature, NOMAD, and scale-factor search
controls. `npcdens()` should keep its argument list focused on estimation
inputs and clearly state that bandwidth-selection controls are accepted
through `...` when bandwidths are omitted and computed internally.

The `.Rd` pilot must not reorder R function formals, alter `\usage{}`
signatures, change defaults, or change examples semantically. Any
reorganization is limited to rendered documentation. Argument group headings
must use the Rd-valid empty sibling form `\subsection{Group Name}{}` followed by
top-level `\item{arg}{...}` entries.

## Pilot Grouping

For `npcdensbw()`:

- main data and bandwidth inputs;
- bandwidth criterion and representation;
- continuous kernel and support controls;
- numerical search controls;
- bounded `cv.ls` quadrature controls;
- categorical kernel controls;
- local-polynomial and NOMAD controls;
- formula-interface controls;
- additional pass-through arguments.

For `npcdens()`:

- main bandwidth and estimation data inputs;
- evaluation and returned-estimate controls;
- pass-through bandwidth-selection controls.

## Gallery Link Policy

The package vignettes point to:

- `https://jeffreyracine.github.io/gallery/np_npRmpi.html`
- `https://jeffreyracine.github.io/gallery/quickstarts.html`
- `https://jeffreyracine.github.io/gallery/kernel_primer.html`
- `https://jeffreyracine.github.io/gallery/density_distribution_quantiles.html`

This pilot may add these as supplemental learning links only after the core
help page remains self-contained.
