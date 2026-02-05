# Issue #3: Local linear with higher‑order kernels shows jumps

**Status:** Open / likely inherent to higher‑order kernel LL

## Repro
```r
library(np)
set.seed(12345)
n <- 1000
x <- runif(n, min=-1, max=1)
y <- x^3 + rnorm(n)

g <- npreg(y~x, regtype="ll", ckerorder=4)
plot(x, fitted(g))
```

## Findings
- Jumps are in **fitted values** (LL intercept), not gradients.
- Higher‑order kernels (4th/6th/8th) allow negative weights.
- LL moment matrix can be **indefinite** even when invertible → discontinuities.
- Several numerical safeguards (determinant checks, ridge, LC fallback triggers) reduced spikes but did not fully remove discontinuities without materially changing estimator.

## Conclusion
Likely a **feature** of higher‑order kernels + LL in finite samples rather than a code defect. A fully smooth fit would require estimator changes (e.g., continuous blending toward LC or boundary rules), which may not be desirable without methodological justification.

## Suggested GitHub response (draft)
- Acknowledge reproducibility.
- Explain higher‑order kernels have negative lobes; LL can become numerically unstable.
- Offer workarounds: use 2nd‑order kernels, local‑constant, or increase bandwidth / apply boundary correction.
