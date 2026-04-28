# `npplreg` Rd Mental Model Audit (2026-04-28)

Scope: `man/np.plregression.Rd` and `man/np.plregression.bw.Rd`.

The partially linear pages now use real `\arguments{}` subheaders with ordinary
top-level `\item{}` entries. The intended reader path is:

1. Identify the linear `X` data, response, nonparametric `Z` data, and
   bandwidth inputs.
2. Choose local-polynomial/NOMAD controls for the nonparametric component only
   when those controls are needed.
3. Adjust search and feasibility controls for component bandwidth searches.
4. Use formula-interface arguments only through the formula method.
5. Pass kernel/support controls through `...` to the component `npregbw`
   searches rather than looking for duplicated formal arguments on
   `npplregbw`.

The estimator page now explicitly tells users that, when `bws` is omitted,
`npplreg()` computes bandwidths through `npplregbw()` and forwards bandwidth
selection controls through `...`. No function signatures or defaults were
changed.
