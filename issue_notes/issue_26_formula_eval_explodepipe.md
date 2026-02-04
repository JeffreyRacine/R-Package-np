# Issue #26: npplregbw formula objects / dynamic calls

Type: Bug

## Report
Using a stored formula or dynamically-constructed formula fails in
`npplregbw()` (and `npcdist()`), due to `explodePipe()` using the raw call
from `match.call()` without evaluating it.

## Reproduction
```
library(np)
df <- data.frame(y=rnorm(10), x1=rnorm(10), x2=rnorm(10), x3=rnorm(10), x4=rnorm(10), x5=rnorm(10))
npplregbw(formula = as.formula(paste0(paste("y ~ ", paste0('x', 1:4, collapse= " + ")), '|x5')), data = df)

fml <- y ~ x1 + x2 + x3 + x4 | x5
npplregbw(formula = fml, data = df)
```

## Proposed Fix
Evaluate non-formula inputs in `explodePipe()` before `as.character()`:
```
if (!inherits(formula, "formula"))
  formula <- eval(formula, parent.frame())
```

## Verification
Both dynamic and stored formula cases run after this change.

## Risk
Low. Only affects formula parsing; evaluated in caller frame to preserve
intended behavior.

## Status
Resolved (evaluation of non-formula inputs before parsing).
