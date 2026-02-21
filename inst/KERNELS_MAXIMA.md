# `kernels.max` Developer Notes

This file is a Maxima script used for symbolic derivation of higher-order
kernel objects (kernel, derivative, CDF, and convolution expressions) and
associated constants.

## Purpose

- `inst/kernels.max` is **developer tooling**, not runtime package code.
- It is useful when extending or validating higher-order kernel algebra.
- It is not sourced by R code during normal package execution.

## Typical Use

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
