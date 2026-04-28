# `npreg` Rd Mental Model Audit (2026-04-28)

Scope: `man/np.regression.Rd` and `man/np.regression.bw.Rd`.

The regression pages now use real `\arguments{}` subheaders with ordinary
top-level `\item{}` entries. The intended reader path is:

1. Identify data and bandwidth inputs.
2. Decide whether the fit is local-constant/local-linear or explicit
   local-polynomial/NOMAD.
3. Choose bandwidth criterion and representation.
4. Adjust search initialization, kernels, and support only when defaults need
   to be changed.
5. Use formula-interface arguments only through the formula method.

The estimator page now explicitly tells users that, when `bws` is omitted,
`npreg()` computes bandwidths through `npregbw()` and forwards bandwidth
selection controls through `...`. This keeps the estimator page concise while
making the bandwidth selector page the canonical guide to the full search
surface.

The bandwidth-selector page separates local-polynomial/NOMAD controls from the
shared search/kernel/support group. All continuous and categorical kernel
arguments remain under `Search Initialization, Kernels, And Support`. No
function signatures or defaults were changed.
