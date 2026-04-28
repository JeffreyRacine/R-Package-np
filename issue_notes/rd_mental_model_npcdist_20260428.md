# npcdist Rd Mental Model Pilot For npRmpi

Date: 2026-04-28

Scope: documentation-only cleanup for `npcdist()` and `npcdistbw()` in
`npRmpi`, following the validated `npcdens` argument-subheader pattern while
preserving MPI lifecycle wording.

## User Decision Sequence

1. Initialize and manage MPI execution as appropriate for the session.
2. Choose formula, data-frame, or precomputed bandwidth interface.
3. Decide whether `npcdist()` should use supplied `bws` or compute bandwidths
   internally through `npcdistbw()`.
4. Choose bandwidth criterion and representation.
5. Choose kernel and support controls, including continuous `cxker*`/`cyker*`
   and categorical `uxkertype`, `oxkertype`, and `oykertype`.
6. Tune numerical search and scale-factor controls only when needed.
7. Use local-polynomial and NOMAD controls only when the polynomial-adaptive
   route is intended.

## Documentation Contract

Use Rd-valid argument subheaders of the form `\subsection{Group Name}{Short
sentence.}` followed by top-level `\item{arg}{...}` entries. Do not change R
formals, `\usage{}`, defaults, examples semantically, or MPI execution
semantics.
