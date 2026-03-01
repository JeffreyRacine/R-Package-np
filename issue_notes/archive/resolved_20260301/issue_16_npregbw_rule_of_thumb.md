# Issue #16: Rule-of-thumb for npregbw

**Status:** Guidance only (no code change)

## Notes
General multivariate mixed‑type regression has no universally valid plug‑in / ROT bandwidth. Users can supply manual bandwidths or use existing CV methods; univariate plug‑ins (e.g., KernSmooth::dpill) are not generally applicable.
