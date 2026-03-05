# `kernels.max` Developer Notes

This file is a Maxima script used for symbolic derivation of higher-order
kernel objects (kernel, derivative, CDF, and convolution expressions) and
associated constants.

## Purpose

- `inst/kernels.max` is **developer tooling**, not runtime package code.
- It is useful when extending or validating higher-order kernel algebra.
- It is not sourced by R code during normal package execution.

## Typical Use

Install/launch Maxima from the official project page:
<https://maxima.sourceforge.io/>

1. Open Maxima.
2. Load `kernels.max`.
3. Run `computeKernels(<base_kernel>, <support>)` for the desired kernel.
4. Capture symbolic/numeric output for:
   - order-2/4/6/8 kernel forms,
   - derivatives,
   - CDF forms,
   - convolution forms,
   - kernel constants.

Example Maxima session sketch:

```maxima
load("kernels.max")$
computeKernels(epanechnikov, sqrt(5))$
computeKernels(gaussian, inf)$
```

## Maintenance Guidance

- Treat outputs as algebraic reference material.
- Before integrating derived expressions into C/R code, verify:
  - numerical parity against existing implementations,
  - boundary/support behavior,
  - objective-value parity in benchmark smoke tests.
- Keep any production code changes separate from symbolic scratch edits.

## Scope Boundaries

- Keep this file for reproducibility of derivations.
- Do not rely on it as an executable part of package installation/runtime.

## Canonical Implementation Directive (2026-03-05)

This repository follows a strict canonical execution rule:

1. One canonical implementation per method (outside explicit `np.tree` branching).
2. Unsupported configurations must fail fast with explicit `stop(...)` diagnostics.
3. No silent remap/coercion of user-selected options (for example `bwmethod`, `regtype`, kernels, `cv.iterate`, or bounds transforms).
4. No hidden alternate execution paths for the same method semantics.
5. All fit-defining options (for example `degree`, `basis`, `bernstein.basis`, kernels, and bounds) must be propagated and used by the canonical path.
6. `np.tree=FALSE` is the default; when `np.tree=TRUE`, behavior must remain semantics-preserving and option-compatible with the canonical path.
7. Remove or reject legacy/debug compatibility branches that add redundant runtime overhead once canonical behavior exists.
