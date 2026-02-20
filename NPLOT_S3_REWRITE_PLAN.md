# npplot S3 Rewrite Plan

## Goals

1. Preserve user-facing behavior for plotting across `np` estimators.
2. Replace monolithic `npplot.*` implementations with class-specific S3 `plot.*` methods plus shared internal helpers.
3. Keep `np-master` plotting core MPI-agnostic so the same core can be reused in `np-npRmpi` with only thin MPI wrappers at method entry.
4. Make future option additions easier by centralizing argument normalization and rendering primitives.

## Current Architecture (Inventory)

### Dispatch chain in current code

- Public call: `plot(<np object>, ...)`
- Current wrappers route into `npplot(bws = x$bws, ...)` (or `x$bw` for partially linear).
- `npplot()` is a generic:
  - `npplot.rbandwidth`
  - `npplot.scbandwidth`
  - `npplot.plbandwidth`
  - `npplot.bandwidth`
  - `npplot.dbandwidth`
  - `npplot.conbandwidth`
  - `npplot.condbandwidth`
  - `npplot.sibandwidth`

### Primary non-`npplot` plot special case

- `plot.npdensity` has a direct 1D fast/display branch (`trainiseval && ndim==1`) before falling back to `npplot`.

### Function size snapshot (monolithic complexity)

- `npplot.rbandwidth`: ~695 lines
- `npplot.scbandwidth`: ~822 lines
- `npplot.plbandwidth`: ~775 lines
- `npplot.bandwidth`: ~645 lines
- `npplot.dbandwidth`: ~630 lines
- `npplot.conbandwidth`: ~909 lines
- `npplot.condbandwidth`: ~912 lines
- `npplot.sibandwidth`: ~496 lines

## Capability Matrix (What Must Be Preserved)

### Cross-cutting options currently supported

- Output behavior: `plot.behavior = "plot" | "plot-data" | "data"`
- Error construction:
  - `plot.errors.method = "none" | "bootstrap" | "asymptotic"`
  - `plot.errors.center = "estimate" | "bias-corrected"`
  - `plot.errors.type = "pmzsd" | "pointwise" | "bonferroni" | "simultaneous" | "all"`
  - `plot.errors.alpha` validation: must lie in `(0, 0.5)`
  - `plot.errors.style = "band" | "bar"`
  - `plot.errors.bar = "|" | "I"`
  - `plot.errors.bar.num`
- Grid / trim / quantile controls:
  - `neval`, `xq`, `yq`, `zq`, `xtrim`, `ytrim`, `ztrim`
- 3D display controls:
  - `perspective` (`persp` path), `theta`, `phi`, `view = "rotate" | "fixed"`, `ticktype`
- Plot styling:
  - `main`, `xlab`, `ylab`, `zlab`, `col`, `lty`, `lwd`, `cex.*`, `xlim`, `ylim`, `zlim`, `type`, `border`, `sub`
- Scaling/layout:
  - `common.scale`
  - global option `plot.par.mfrow`

### Estimator-specific capabilities

- `rbandwidth` (regression): gradients support; GLP derivative order support (`gradient.order`).
- `scbandwidth` (smooth coefficient): coefficient-function plotting and optional gradients.
- `plbandwidth` (partially linear): component plotting including nonparametric part and scale by parametric coefficients.
- `bandwidth` / `dbandwidth` (unconditional density/distribution families): density/distribution plotting with errors.
- `conbandwidth` (conditional density): conditional surfaces/slices and gradients.
- `condbandwidth` (conditional distribution): conditional CDF surfaces/slices, quantile-regression mode flag handling.
- `sibandwidth` (single-index): single-index curves, gradients, bootstrap/asymptotic errors.

## Duplication to Modularize

## High-duplication patterns in `np.plot.R`

1. Repeated argument normalization/validation blocks:
   - `match.arg` chains
   - `plot.errors.alpha` range checks
   - bootstrap method/block-length defaulting
2. Repeated error-range construction:
   - `compute.default.error.range`
   - `compute.all.error.range`
   - repeated in-method switch logic for `plot.errors.type == "all"`
3. Repeated draw logic:
   - core line plotting + optional bands/bars
   - repeated legend creation for `"all"` bands
4. Repeated data-return assembly when `plot.behavior != "plot"`
5. Repeated 2D/3D branching and fixed/rotate view behavior.

## Proposed Internal Modules

1. `.np_plot_normalize_args()`
2. `.np_plot_validate_errors()`
3. `.np_plot_build_eval_grid()`
4. `.np_plot_compute_errors()` (method-specific adapter to bootstrap/asymptotic engines)
5. `.np_plot_render_2d()`
6. `.np_plot_render_3d()`
7. `.np_plot_render_errors()`
8. `.np_plot_collect_output()` (for `plot.behavior`)
9. `.np_plot_legend_all_bands()`

## Target S3 Design

## Public layer

- Use class-specific `plot.*` methods as first-class entry points.
- `npplot` remains temporarily as compatibility shim only.
- `plot.*` methods call estimator-specific internal engines, not `npplot.*` monoliths.

## Internal layer

- Thin estimator engines:
  - `.np_plot_npregression()`
  - `.np_plot_npcdensity()`
  - `.np_plot_npcdist()`
  - `.np_plot_npudensity()`
  - `.np_plot_npudist()`
  - `.np_plot_npscoef()`
  - `.np_plot_npplreg()` / `.np_plot_singleindex()` as needed
- Shared modules (above) handle common normalization, error logic, and rendering.

## Documentation Plan

1. Keep `?npplot` during transition but mark as internal compatibility gateway.
2. Move user-facing plot argument documentation to class-relevant `plot.<class>` docs.
3. Keep one shared doc topic for common plotting argument semantics to prevent drift.
4. Ensure examples remain `plot(<fit>)` focused; avoid direct `npplot()` examples in user docs.

## Migration Plan (Phased)

### Phase 0: Freeze compatibility contract

- Build a parity checklist covering:
  - estimator class
  - 2D/3D mode
  - error method/type combos including `"all"`
  - gradients and GLP derivative order behavior
  - returned object shape for `plot.behavior`

### Phase 1: Extract cross-cutting helpers

- Introduce internal helper modules with no behavior change.
- Switch existing `npplot.*` methods to call helpers.

### Phase 2: Introduce estimator-specific `plot.*` engines

- Rewire `plot.<np object class>` wrappers to call new engines directly.
- Keep `npplot.*` methods as wrappers to new engines for backward compatibility.

### Phase 3: Retire `npplot.*` monoliths

- Remove old bodies once parity tests are green.
- Keep minimal `npplot` compatibility if needed for one release cycle.

### Phase 4: Port to `np-npRmpi`

- Port shared plotting core as-is.
- Add MPI-specific guarding/dispatch only at public `plot.*` entry points.
- No MPI logic inside rendering helpers.

## Test/Validation Requirements

1. Snapshot path coverage (at least smoke-level) for each estimator family:
   - `perspective = TRUE/FALSE`
   - `plot.errors.method = none/bootstrap/asymptotic`
   - `plot.errors.type = pmzsd/pointwise/bonferroni/simultaneous/all`
   - `plot.behavior = plot/plot-data/data`
2. GLP checks:
   - `regtype="glp"` with derivative orders and consistency against existing behavior.
3. Error validation checks:
   - `plot.errors.alpha` rejects values outside `(0,0.5)` with current error text.
4. `npRmpi` checks:
   - no-hang smoke for supported plot paths.

## Explicit Non-Goals (for this rewrite)

1. No `ggplot2` migration.
2. No changes to estimator math.
3. No API break to existing `plot()` calls.

## Risks and Mitigations

1. Risk: subtle output regressions in `"all"` error-band rendering.
   - Mitigation: dedicated parity tests for all-band overlays and returned data objects.
2. Risk: duplicated semantics drift between docs and code.
   - Mitigation: centralize common plot argument docs and reuse references.
3. Risk: accidental MPI coupling in core plotting helpers.
   - Mitigation: enforce MPI wrapper-at-entry-only rule during review.
