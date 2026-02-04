# Issue #7: uocquantile() with factor subsets

**Type:** Bug

## Report
`uocquantile()` returns an incorrect mode when called on a factor subset that drops the original mode. Example from issue:

```
fruit.l <- as.factor(c("a","n","a","n","a","s"))
uocquantile(fruit.l,.5)              # a
uocquantile(fruit.l[fruit.l!="a"],.5) # expected n, returns s
```

Cause: `table(x)` includes unused levels, while `unique(x)` does not, so the index `j` can map to a wrong level.

## Reproduction (local)
- With baseline code: `uocquantile(fruit.l[fruit.l!="a"], .5)` returns `s`.
- After patch: returns `n` (correct).

## Fix
In `R/np.plot.R`, `uocquantile()` now uses `droplevels(x)` for factors/ordered factors and maps `j` to `levels(x)`. This ensures the counts and labels align.

### Patch summary
- `x <- droplevels(x)` for factor/ordered.
- Use `levels(x)` (instead of `sort(unique(x))`) to map the mode index.

## Risk assessment
Low risk. Only affects factor/ordered paths and corrects behavior when unused levels are present. For non-factor data, unchanged.

## Status
Implemented locally in `R/np.plot.R`. Needs package install to propagate.

## Status
Resolved (factor levels dropped before mode computation).
