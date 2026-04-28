# npcdist Rd Mental Model Pilot

Date: 2026-04-28

Scope: documentation-only cleanup for `npcdist()` and `npcdistbw()` in `np`,
following the validated `npcdens` argument-subheader pattern.

## User Decision Sequence

1. Choose formula, data-frame, or precomputed bandwidth interface.
2. Decide whether `npcdist()` should use supplied `bws` or compute bandwidths
   internally through `npcdistbw()`.
3. Choose bandwidth criterion and representation.
4. Choose kernel and support controls, including continuous `cxker*`/`cyker*`
   and categorical `uxkertype`, `oxkertype`, and `oykertype`.
5. Tune numerical search and scale-factor controls only when needed.
6. Use local-polynomial and NOMAD controls only when the polynomial-adaptive
   route is intended.

## Documentation Contract

Use Rd-valid argument subheaders of the form `\subsection{Group Name}{Short
sentence.}` followed by top-level `\item{arg}{...}` entries. Do not change R
formals, `\usage{}`, defaults, or examples semantically.
