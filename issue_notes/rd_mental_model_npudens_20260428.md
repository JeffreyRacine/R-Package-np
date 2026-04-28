# `npudens` Rd Mental Model Audit (2026-04-28)

Scope: `man/np.density.Rd` and `man/np.density.bw.Rd`.

The unconditional density pages now use real `\arguments{}` subheaders with
ordinary top-level `\item{}` entries. The intended reader path is:

1. Identify data and bandwidth inputs.
2. Choose the bandwidth criterion and representation.
3. Adjust search initialization, kernels, and support only when defaults need
   to be changed.
4. Use formula-interface arguments only through the formula method.

The estimator page now explicitly tells users that, when `bws` is omitted,
`npudens()` computes bandwidths through `npudensbw()` and forwards bandwidth
selection controls through `...`. This keeps the estimator page concise while
making the bandwidth selector page the canonical guide to the full search
surface.

The bandwidth-selector page keeps continuous kernel, categorical kernel,
support, and search-start controls in one group because they jointly determine
the feasible bandwidth-search problem. No function signatures or defaults were
changed.
