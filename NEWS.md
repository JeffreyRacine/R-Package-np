# npRmpi 0.70-6

* Conditioned local-polynomial influence and response owners now select the
  same globally stable graded basis before applying the canonical diagonal
  ridge policy. The former influence-only empirical QR rotation spanned the
  same polynomial space, but diagonal ridging is not invariant to that extra
  change of representation; on translated or calendar-scale predictors it
  could therefore make `npreghat()` matrix/apply disagree with `npreg()` even
  though both used the same factor policy. Objective, fit, matrix, apply, and
  derivative consumers now share the fit-capable representation selected by
  one helper. Positive-degree `npreghat()` and `npcdhat` matrix/apply routes
  now enter that same regression core; the retired QR owner, private
  positive-degree R kernel-weight owner, conditional per-row owner, and their
  workspaces have been removed. Raw well-conditioned arithmetic is unchanged,
  no collective was added, and auxiliary storage remains linear in the sample
  size.

* Canonical wider-than-one local-polynomial ridge retries are now invariant to
  a common scaling of all kernel weights. The shared response/adjoint factor
  owner interprets `1/n.train` as a dimensionless fraction of the pristine
  Gram scale, using one allocation-free diagonal scale for every bounded
  retry. Fixed-objective, fit, leave-one-out apply, and hat owners therefore
  select the same regularized estimator even when one kernel provider retains
  the bandwidth normalization and another omits the cancelling common factor.
  Zero-ridge arithmetic is unchanged; formerly ill-conditioned objective and
  hat values can change to the scale-equivariant canonical result.

* Fixed local-polynomial CVLS, check-loss, and Klein--Spady objectives no
  longer divert rows with at most one donor per basis term through a private
  `DGELSY` minimum-norm fit. Those rows now use the same condition-qualified,
  bounded canonical response solve as every other fixed row and as public
  leave-one-out matrix/apply. Removing the retired fork also removes its
  observation-by-basis support-index, weight, and count workspaces, pair-loop
  branches, and MPI support-only reductions. Objective values can change for
  underdetermined compact-support rows that formerly entered the distinct
  minimum-norm estimator.

* Fixed-resident and adaptive-BLAS local-polynomial regression objectives now
  use the same condition-aware response/factor policy as ordinary fits,
  leave-one-out influence rows, and the independent hat-matrix owner. CVAIC
  obtains leverage from the retained factor through the canonical adjoint
  solve, so response and influence calculations share one bounded ridge
  transcript instead of selecting regularization independently. The MPI
  fixed-LP fit owner follows the same response policy without adding a
  collective. This removes the remaining LU-failure-only ridge loops;
  well-conditioned CVLS arithmetic is unchanged, while formerly
  ill-conditioned objective and fit results can change to the canonical
  condition-qualified solution.

* Positive-degree local-polynomial fits now share the conditional objective's
  canonical, basis-neutral coordinate policy across coordinator and rank-owned
  routes. Raw polynomial coordinates are retained when well conditioned;
  otherwise ordinary and all-large fit owners select one stable global
  representation for training, evaluation, derivatives, and uncertainty
  calculations. The shared design validator now assesses intrinsic rank in
  those same stable coordinates, so calendar-scale predictors are no longer
  mistaken for rank-deficient designs while genuine collinearity still fails.
  This is a representation change within the same polynomial span, not
  evaluation-point centering or an unbounded ridge fallback, and uses only
  linear-in-sample basis storage without adding MPI communication.

* Rank-local generalized nearest-neighbour regression evaluations now share
  one execution-context classifier across `npreg()` and `npreghat()`. This
  prevents a worker from re-entering a nested MPI kernel-weight owner merely
  because the active-pool flag is coordinator-local, preserving persistent
  pools across successive kernel and local-polynomial calls without adding a
  collective, changing numerical arithmetic, or introducing size-based
  dispatch.

* Ordinary univariate generalized nearest-neighbour regression now admits
  `k = 1` for the positive second-order Gaussian and Epanechnikov LC mean and
  derivative operators, including the equivalent generalized-LP degree-zero
  representation. One R-owned capability resolver supplies manual validation,
  optimizer bounds, and the native lower bound; adaptive NN, mixed or
  multivariate regressors, positive LP degree, higher-order/uniform/bounded
  kernels, and semiparametric owners retain `k >= 2`. Literal zero radii remain
  infeasible objectives and explicit fit/hat errors through the shared
  occurrence-aware radius validator. Fixed-degree searches evaluate the newly
  admitted discrete endpoint once, while automatic degree searches whose grid
  includes a positive degree retain the safe joint `k >= 2` capability envelope.
  The resolver is regression-scoped; density and distribution bandwidth
  classes retain their existing nearest-neighbour domains.

* Smooth-coefficient NOMAD degree search now transports its scoped prepared
  context through an exact locked opcode token instead of comparing a
  serialized function closure by process-local identity. Its service workers
  are also stopped before the context lease is retired. This restores
  one- and multi-worker `npscoefbw()` completion without changing objective
  arithmetic, search policy, evaluation accounting, or message count.
* Regression leave-one-out hat-matrix and apply routes now share the objective's
  occurrence-aware generalized/adaptive NN geometry, source basis, weighted
  moment layout, factorization, bounded ridge transcript, and intercept
  restoration. Matrix and apply consume the same completed influence row;
  apply retains linear-in-sample auxiliary storage and does not construct an
  evaluation-by-training matrix. The canonical factor owner now treats an LU
  factor below the working-precision reciprocal-condition threshold as
  numerically rank deficient before applying the existing bounded ridge rule.
  Well-conditioned `k >= 2` objectives remain pointwise unchanged; formerly
  unstable rank-deficient results can change, and all-zero compact-kernel LOO
  rows now fail explicitly instead of returning a fabricated ridged estimate.

* Exact adaptive-nearest-neighbour leave-one-out objectives now construct each
  donor radius on the literal fold that excludes both the focal observation
  and the donor itself. A shared adjacent-order-statistic owner supplies
  regression, unconditional density, conditional density, and empirical
  conditional-distribution rows without a per-fold sort or pair matrix.
  Saturated and extended counts scale the largest available fold radius by
  `k/(n - 2)`; the full-design integrated-square term in density CVLS retains
  its distinct full-sample radius and `k/(n - 1)` extension. This repairs
  formerly invalid high-count objectives and adaptive regression optimizer
  handoffs while preserving literal duplicate occurrences, zero-radius
  infeasibility, and linear-in-sample auxiliary storage.

* Extended generalized-nearest-neighbour training queries now preserve the
  supplied occurrence identity before applying the package's linear radius
  extension. Regression, unconditional density/distribution, and conditional
  mapped-training owners therefore select the farthest literal neighbour only
  after deleting the identified focal occurrence. Equal-valued external
  queries retain external geometry, and ordinary non-extended counts are
  pointwise unchanged.

* Empirical-sample conditional-distribution CV now evaluates the documented
  off-diagonal criterion with divisor `n(n - 1)`. Ordinary generalized
  nearest-neighbour rows construct the explanatory radius from the delete-one
  sample and the response radius from the occurrence union `{i, j}` required
  by `Fhat_{-i}(Y_j | X_i)`. The exact row and bounded block owners share one
  two-slot order-statistic selector, preserve duplicate occurrences, and add
  no pair matrix or per-fold sort. Apart from the separately documented exact
  adaptive-fold correction, fixed, adaptive, and beta sample-grid objectives
  receive only the normalization correction; genuine external and default
  response grids are unchanged. Because the diagonal was already
  omitted, this changes reported training-grid criterion values by the common
  positive factor `n/(n - 1)` relative to the former `n^2` divisor and does
  not change the mathematical minimizer on a common candidate domain.

* Empirical distribution and conditional-distribution CVLS owners now
  finalize their shared pair-count normalization through one checked native
  seam. Invalid training/evaluation shapes and non-finite accumulators return
  the existing infeasible status instead of permitting a zero pair-count
  division to produce `NaN` or `Inf`; valid objective arithmetic, MPI
  reductions, search domains, and public interfaces are unchanged.

* Ordinary generalized nearest-neighbour conditional-density and
  conditional-distribution training fits now construct explanatory and
  response radii after deleting the identified focal occurrence, while still
  retaining that observation in the estimator sum. Scalar and
  local-polynomial owners share the same explicit occurrence contract and
  reject literal zero radii. Equal-valued external queries, fixed and adaptive
  bandwidths, beta kernels, extended nearest-neighbour bandwidths, and
  external-grid search objectives retain their existing contracts.

* Ordinary generalized nearest-neighbour conditional-density CVML now
  constructs both explanatory and response radii after deleting the focal
  training occurrence used by each leave-one-out objective row. The scalar,
  local-polynomial row, and bounded block-stream owners share the same
  explicit occurrence contract; zero literal radii fail the candidate
  cleanly. Fixed and adaptive bandwidths, conditional-density CVLS,
  conditional-distribution external-grid objectives, beta kernels, and
  extended nearest-neighbour bandwidths retain their existing contracts.

* Ordinary generalized nearest-neighbour unconditional density and
  distribution rows now distinguish identified training occurrences from
  equal-valued external evaluation points. Training fits and CVML construct
  each query radius after deleting its focal occurrence, while full fitted
  rows continue to use every observation in the kernel sum. External queries,
  density CVLS, adaptive and extended nearest-neighbour bandwidths, and beta
  kernels retain their existing contracts; a literal zero radius is rejected
  rather than replaced by a positive neighbour.

* Broad-support positive-degree `npscoef` tree objectives now retain one
  rank-local, invocation-owned outer-product workspace across full-contiguous
  tree rows. The unchanged design matrix is packed once while response/weight
  scratch is reused; sparse, partial, and multi-segment support retains the
  incumbent row-local owner. Public bandwidths, objectives, evaluation
  accounting, fitted results, and coarse MPI row ownership are unchanged.

* `nplsqregbw()` and `nplsqreg()` now support the package-wide
  `nomad = "auto"` shortcut. One continuous smoothing dimension uses a
  family-local coordinator that gives each exhaustive fixed-degree candidate
  its own complete autodispatch transaction; two or more dimensions retain
  the existing NOMAD/Powell joint search. Historical `nomad = FALSE` and
  `nomad = TRUE` routes, the shared MPI objective owner, and the separate
  logical-only `nomad.pilot` control are unchanged.

* Conditional-density automatic degree search now provisions response-side
  categorical workspace for every engine that a later degree candidate can
  activate, rather than only the engine selected by the initial degree. This
  restores `npcdensbw()` and `npconmode()` search-trajectory parity across
  serial and MPI owners without changing fixed-candidate objective arithmetic,
  degree policy, or public evaluation accounting.

* Regression bandwidth plots now request standard errors on every direct
  asymptotic-fit call used by the MPI-aware progress path. This preserves the
  shared plotting contract for unavailable gradient panels, endpoint geometry,
  and ordinary formula/session plots instead of failing during interval
  assembly.

* On Apple silicon, fixed-bandwidth compact-support regression-tree CV
  objectives of local-polynomial basis width six now use compile-time
  resident sparse-pair and rank-owned-row microkernels. Tree order, coarse MPI
  outer-row ownership, the unique upper Gram triangle, and the common
  mirror/solve finisher are unchanged. Widths below and above six are
  unaffected, and portable builds retain the incumbent generic width-six
  owner.

* Fixed-bandwidth compact-support regression CV candidates whose kernel
  support covers the complete observed continuous range now use the canonical
  dense LP owner instead of unproductive full-tree traversal. A shared,
  operator-aware geometry classifier also routes full-support fixed kernel
  sums to their dense sibling and searches mixed-support trees only over
  coordinates whose kernels can prune, while continuing to evaluate every
  kernel factor. These exact, allocation-free, rank-symmetric policies apply
  to Apple and portable builds without sample-size or runtime thresholds and
  without new communication; width-two regression CV, genuinely sparse,
  nearest-neighbour, non-compact, and unsupported owners retain their
  established routes.

* Conditional-density native degree search now constructs its finite invalid-
  candidate penalty once from the same canonical baseline-and-retry rule as
  the ordinary bandwidth route and restores that rank-common value before
  every candidate. This repairs false `-10000000` NOMAD baselines in
  `npcdens()` and `npconmode()` while preserving the incumbent vectorized
  width-one CVML owner, public evaluation accounting, and linear rank-local
  storage.

* MPI autodispatch now materializes S3 forwarded arguments exclusively from
  the method frame that owns their promises. Named `...` controls, explicit
  `NULL`, and omitted arguments retain their identities without scanning
  unrelated dynamic frames. This repairs stack-dependent NOMAD multistart
  transport in `nplsqregbw()` and removes its obsolete route-local frame
  assignment while preserving the documented public defaults and native
  search semantics.

* `options(np.macMseries.accelerate = FALSE)` now disables the
  Apple-qualified weighted-design and conditional LP transpose-GEMV siblings
  consistently. Ordinary BLAS/LAPACK algebra remains available independently
  of this option, and `TRUE`/`"auto"` retain the qualified Apple-arm64 owners.

* Generalized fixed-bandwidth local-polynomial CV objectives of basis width
  six now retain only the unique upper triangle of their symmetric weighted
  Gram rows and mirror it once before solving. Serial unordered-pair and MPI
  rank-owned rows share this arithmetic contract; the portable scalar and
  Apple-arm64 vector siblings avoid the packed-width regression when optional
  acceleration is disabled while preserving wider, tree, and nearest-neighbour
  owners and linear-in-sample auxiliary storage.

* `npconmode()` now evaluates categorical response levels through one
  memory-bounded conditional-density batch instead of paying a complete MPI
  fit lifecycle for every level. Requested class-probability gradients retain
  one singleton evaluation for their selected level while all other levels
  share the batch. The canonical rank-owned `npcdens()` engine, probability
  repair, class selection, and result arithmetic are unchanged; temporary
  evaluation rows are capped independently of the number of response levels.

* Eligible fixed-bandwidth local-polynomial CV objectives now apply the same
  basis-neutral owner rule to generalized, additive, and tensor bases:
  resident accumulation through basis width six, packed BLAS from width seven,
  and the independently qualified generic resident owner for sparse tree rows.
  The hoisted basis, canonical solve, deletion, projection, loss arithmetic,
  and MPI row ownership are shared; no basis-family or kernel-spelling
  fallback remains. Widths five and six retain only the unique symmetric Gram
  triangle during accumulation and mirror it once before solving.

* Fixed-bandwidth compact-support `npscoef` objectives now pack each rank-local
  sparse tree row's active support into the existing BLAS outer-product engine
  for moderate and wide local-polynomial moment systems. Coarse outer-row MPI
  ownership is unchanged, while dense, degree-zero, nearest-neighbour, signed,
  derivative, and permutation routes retain their existing contracts and
  rank-local workspace bounded by active support and moment-system width.

* Fixed-bandwidth compact-support regression trees now use a bounded generic
  resident-row accumulator for local-polynomial basis widths above five,
  sharing the canonical solve, deletion, and uncentered projection engine.
  The MPI owner retains coarse outer-row distribution, so explicit
  `np.tree = TRUE` and eligible `np.tree = "auto"` both use the repaired tree;
  dense and smaller-width routes are unchanged.

* Fixed-bandwidth degree-zero conditional-density CVLS now recognizes an
  all-categorical explanatory side at its constant upper kernel limit as the
  canonical global-X design. This restores the unconditional-density limit,
  including the ordered Li-Racine case, and avoids redundant conditional X
  weighting. Generalized/adaptive NN and response-side kernel limits are
  unchanged.

* Degree-zero conditional-density and conditional-distribution CVLS objectives
  now enter the same canonical all-large explanatory-bandwidth context as
  positive-degree local polynomials. The width-one case uses its implicit unit
  basis and scalar Gram inverse, avoiding the dense locally weighted pass while
  preserving delete-one objective semantics and linear rank-local storage.
  Conditional-distribution tree rows now translate their independent response
  permutation once at the shared indicator boundary.

* Fixed-bandwidth two-predictor local-linear CV objectives now keep sparse
  tree rows resident and accumulate the same unique symmetric Gram triangle as
  the canonical dense width-three engine before using their shared mirror,
  solve, and uncentered projection. Serial pair enumeration and rank-owned MPI
  rows share this algebra; MPI still distributes complete outer rows and
  performs one bounded collective reduction.

* Bootstrap-intensive tests now reuse fixed kernel geometry in bounded,
  rank-local replication groups. `npcmstest()` and `npqcmstest()` keep one
  model refit per draw but contract the resulting residual or score columns
  through multi-response kernel sums; fixed, unbounded `npdeneqtest()` routes
  process pooled multiplicity columns while bounded and nearest-neighbour
  routes keep literal duplicate resampling; and
  `npdeptest(method = "summation")` uses a bounded native index batch that
  shares only the invariant response-side marginal work. Draw order, complete
  bootstrap payloads, bandwidth and model semantics, RNG restoration,
  progress, rank-count independence, and bounded rank-local workspace are
  preserved.

* Spawned MPI sessions now finalize the master during an early process-exit
  finalizer, after closing and verifying any retained worker pool, instead of
  deferring finalization until namespace teardown. Normal scripts consequently
  return status zero without external `setsid` containment, including soft
  pool reuse, forced close, and omitted explicit close. Public pool-reuse and
  `mpi.finalize()` contracts are unchanged; standard tests and demo launchers
  now require a clean exit and no longer accept status 137 or broadly kill
  processes by name.

* Fixed-bandwidth fourth- and sixth-order Gaussian convolution rows now hoist
  their invariant normalization scale out of the observation loop. Conditional
  CVLS objective values are unchanged while the affected higher-order kernel
  pass avoids two divisions per row element.

* Gaussian entropy integration now hoists validated bandwidth reciprocals out
  of its training-by-evaluation kernel loop and evaluates its three
  exponentials in bounded vForce blocks on supported Apple silicon. This
  accelerates `npdeptest()` and `npsdeptest(method = "integration")` without
  changing bandwidths, bootstrap decisions, accumulation order, or asymptotic
  memory use; `options(np.macMseries.accelerate = FALSE)` retains the scalar
  route.

* Default fixed-Gaussian entropy summation now uses registered, streamed
  native owners for `npunitest()`, `npsymtest()`, `npdeptest()`, and
  `npsdeptest()` on every rank. Univariate and symmetry bootstraps consume
  bounded chunks of multiplicities while preserving the exact incumbent
  `sample.int()` or `boot::tsboot()` draw plan, RNG state, collective fanout,
  and reassembly; dependence statistics fuse their marginal and joint kernel
  sums without an n-by-n cache. Explicit generalized- and adaptive-nearest-
  neighbour requests remain on literal, duplicate-preserving resample/refit
  routes so sample-owned radii and ordered lag pairs are unchanged.

* General bounded conditional-density CVLS now reuses each response-side
  quadrature kernel tile across a memory-bounded group of rank-owned
  explanatory-row tiles. Ordinary and beta kernels share the same canonical
  engine; MPI ownership, objective arithmetic, and accumulation order are
  unchanged, and storage remains `O(n B)` under the existing conditional-LP
  tile budget.

* Eligible bivariate surface plots now share one base-perspective frame owner,
  one viridis palette mapper, and one neutral-gray three-face grid style across
  the base and `rgl` renderers. The shared renderer is used by the regression,
  density, distribution, conditional, semiparametric, conditional-mode, and
  copula plot families; estimator and plot-data arithmetic are unchanged.

* Internal direct regression evaluators now enter native MPI collectives only
  when the current call is executing symmetrically on every rank.  Master-only
  fit, prediction, and test helpers remain local even when a slave pool is
  active, preventing an unmatched collective wait; public auto-dispatched and
  explicit broadcast regression calls retain their distributed owner and
  estimator arithmetic.

* Prepared local-polynomial conditional-density maximum-likelihood
  cross-validation now assigns each observation row to one MPI rank and
  reduces an order-preserving contribution vector once per objective. This
  removes rank-replicated LP row accumulation while retaining byte-identical
  objective values and bounded workspace. The incumbent nonprepared route and
  the beta kernel's existing distributed row engine remain isolated in their
  original sibling kernels, avoiding nested ownership and adjacent timing
  regressions.

* Fit and evaluation uncertainty now use one request vocabulary: `se = FALSE`
  for `npreg()`, `npscoef()`, `npindex()`, and `nplsqreg()`;
  `se.fit = FALSE` for
  prediction; and `gradients(..., se = FALSE)` for gradient standard errors.
  The former Boolean estimation argument `errors` is intentionally rejected
  with a migration message, while plot methods retain their method-valued
  `errors` control. Explicit `se = TRUE` reproduces the previous uncertainty
  arithmetic; the default avoids unrequested moments, bootstrap work,
  covariance assembly, and native output vectors. `se()`, `vcov()`, and
  gradient-SE extraction fail helpfully when the required state was not
  computed. Bandwidth objectives, MPI ownership, and plot interval semantics
  are unchanged.

* Conditional density and distribution kernel summaries now identify
  explanatory and dependent kernels separately for continuous, unordered, and
  ordered variables. Bandwidth-object `print()` and `summary()` methods share
  one formatter, public output no longer exposes internal `cxker*`/`cyker*`
  argument mnemonics, and invalid categorical-kernel selections name the
  specific `uxkertype`, `uykertype`, `oxkertype`, or `oykertype` argument.

* During bandwidth-objective evaluation, MPI-3 builds now use progress-aware
  nonblocking forms of the existing shared post-kernel collectives. Rank zero
  can therefore maintain the canonical progress heartbeat while waiting for
  workers, without changing collective buffers, operations, ordering, or
  estimator arithmetic. Fit/evaluation routes and MPI-2 builds retain their
  existing blocking collectives.

* Interactive `renderer = "rgl"` surfaces now use a high-DPI-aware widget
  canvas, a higher-resolution uncertainty legend, and an initial camera more
  closely matched to the base perspective. The renderer caps backing-store
  scaling at two device pixels per CSS pixel, retains native rgl axis labels
  and interaction, and preserves all existing `rgl.*` user overrides. Shared
  rgl setup and option routing now have one internal owner across supported
  plot families; estimator, bootstrap, interval, and base-renderer arithmetic
  are unchanged.

* `npindexbw()` now canonicalizes its numeric design matrix once at the public
  boundary and reuses that matrix during internal index and coordinate setup.
  This removes redundant full-matrix copies before objective evaluation while
  preserving estimator arithmetic and conditions on every rank.

* The integration routes for `npunitest()`, `npsymtest()`, `npdeptest()`, and
  `npsdeptest()` now use bounded deterministic corrected-trapezoid quadrature
  for their default Gaussian statistics. Univariate iid bootstraps reuse a
  fixed grid and evaluate bounded chunks of resampling counts; symmetry
  bootstraps reuse reflected density values; and bivariate quadrature calls
  the registered constant-storage Gaussian evaluator in bounded tiles. This
  removes the `cubature` dependency while preserving bootstrap sampling laws,
  fixed bandwidth reuse, result shapes, seed restoration, nondefault kernel
  routes, and existing MPI fanout and reassembly contracts.

* Removed the private, unexported package-local `gsl.bs` bridge and its
  registered native wrappers. Public spline construction remains owned by
  `crs::gsl.bs()`; the independent B-spline primitives used by the canonical
  Bernstein local-polynomial engine remain in place on every rank, with
  estimator arithmetic unchanged.

* The canonical beta row engine now owns its malloc-backed route, categorical,
  derivative, and row workspaces behind one unwind-protected cleanup boundary
  on every rank. Nested derivative helpers borrow that invocation-scoped
  scratch rather than allocating it per call, so user interrupts and R errors
  cannot strand native storage. Successful PDF, CDF, convolution, derivative,
  powered-derivative, and regression arithmetic and rank ownership are
  unchanged.

* The shared regression fit now allocates its operator vector once, after
  engine and route validation. This removes a duplicate general-path
  allocation and the scalar beta route's retained allocation on every rank
  without changing estimator arithmetic.

* Beta nearest-neighbour bandwidth preparation now calls the shared
  continuous-distance owner directly on every rank. It no longer temporarily
  rewrites unrelated package-global scaling and categorical factors, while
  retaining the same generalized- and adaptive-nearest-neighbour arithmetic.

* Nearest-neighbour lookup and distance failures now leave both shared
  bandwidth owners through their normal cleanup exit on every rank. This
  removes retained standard-deviation and distance workspaces after
  recoverable failures while preserving successful fixed, generalized-
  nearest-neighbour, and adaptive-nearest-neighbour arithmetic and
  conditions.

* Removed the unreachable conditional-bootstrap exact-state sidecar and its
  bypassed MPI worker block. Public conditional bootstrap routes retain their
  canonical distributed refit/count and fixed local-constant owners, without
  serial fallback.

* Removed the unregistered and callerless internal `beta_kernelsum` facade.
  Beta bandwidth-mode ownership now resides in the shared bandwidth contract;
  all estimator, kernel-sum, LP, optional-output, and distributed routes
  continue to use the canonical continuous-row engine.

* Higher-order beta convolution now prepares each distinct side/component
  shape once per observation pair and rank rather than rebuilding it for every
  opposing component. The stack-local decomposition retains the canonical
  centered gamma-ratio and signed component-pair accumulation order, changes
  no order-two arithmetic, and adds no dynamic or sample-sized storage.

* Fixed-bandwidth beta CDF rows now prepare each coordinate/component
  concentration once per invocation and rank rather than recomputing it for
  every evaluation-observation pair. The fixed-only sibling shares the
  canonical CDF arithmetic and failure statuses with scalar and nearest-
  neighbour rows, cannot fall back after selection, and adds only coordinate-
  by-component transient storage per rank.

* Beta CDF rows now prepare observation support coordinates and higher-order
  coefficient state once per invocation and rank rather than rebuilding them
  for every evaluation-observation pair. Fixed, generalized-nearest-neighbour,
  and adaptive-nearest-neighbour distribution rows share the same canonical
  `pbeta()` arithmetic and failure contract. The prepared route is selected
  outside row traversal, cannot fall back to the scalar CDF evaluator, and
  uses only beta-coordinate-linear transient storage per rank; density-only
  routes retain their existing prepared-PDF path unchanged.

* Beta regression gradient rows now form the signed-log level and its
  regular/jump target derivative from one canonical component preparation on
  every rank. The level phase still completes before derivative work,
  preserving native failure precedence and exact fitted, gradient, and
  standard-error arithmetic while removing the duplicate beta shape and PDF
  pass. Fixed- and generalized-nearest-neighbour gradient rows additionally
  reuse their prepared evaluation-row components on every rank.
  Adaptive-nearest-neighbour gradient rows reuse invocation-owned observation
  transforms and coefficient state while retaining pair-owned bandwidth
  shapes and normalizers. Each bandwidth topology selects an isolated sibling
  before row traversal, with no inner-loop route branch or scalar fallback.
  Derivative-only fixed/GNN state is allocated lazily on every rank, so
  level-only fits retain their incumbent preparation cost.

* Adaptive-nearest-neighbour beta PDF rows now prepare support-transformed
  observations and higher-order coefficient state once per invocation while
  retaining observation-owned bandwidth shapes and normalizers in the
  pairwise evaluator. This preserves exact kernel-weight, objective, fit, and
  prediction arithmetic on every rank, uses beta-coordinate-linear transient
  storage per rank, and cannot fall back to the scalar row after the prepared
  route is selected.

* Fixed- and generalized-nearest-neighbour beta PDF rows now prepare support-
  transformed observations, component normalizers, and higher-order
  coefficient logs at their invocation or evaluation-row lifetime instead of
  recomputing them for every observation pair. The shared row engine retains
  exact objective, fit, prediction, and kernel-weight arithmetic, has no
  runtime fallback after route selection, and uses only beta-coordinate-linear
  transient storage per rank.

* Hardened beta self-map layout validation at the registered kernel-sum,
  nearest-neighbour bandwidth, and canonical continuous-row boundaries.
  Internal calls that identify training and evaluation data now require equal
  extents and a strict logical self-map flag, preventing an inconsistent
  native invocation from aliasing a shorter training buffer as a longer
  evaluation buffer. Valid estimator and `npksum()` arithmetic is unchanged.

* Corrected automatic-degree regression startup for beta kernels and
  local-polynomial engine transitions for every continuous kernel. Internal
  zero bandwidth placeholders are now owned by bandwidth selection rather
  than validated as manual beta bandwidths, and each evaluated degree selects
  the scalar width-one or general LP engine through the shared canonical rule
  on every rank. Returned objectives replay exactly across one- and
  multi-worker searches without fallback or post-search rewriting.

* Corrected native regression NOMAD handoff so the returned bandwidth vector
  is the exact decoded point evaluated by the resident objective, including
  mixed predictors stored in non-kernel-grouped column order. The resident
  search now shares the ordinary regression kernel-state configuration for
  beta order, categorical compression, and tree eligibility; fixed-degree and
  automatic-degree results therefore replay their stored objective exactly
  across local, one-worker, and multi-worker routes.

* Corrected fixed-bandwidth beta local-polynomial regression with an active
  MPI pool so rank ownership retains the canonical scaled beta moment row
  used by serial and nearest-neighbour routes. The owner path no longer
  substitutes the legacy continuous-kernel accumulator; fitted values,
  standard errors, and gradients now agree across local, one-worker, and
  multi-worker evaluation while retaining bounded row storage.

* `npseed()` now applies its validated C-backend seed to every rank of an
  active MPI pool through the canonical command-broadcast path. This corrects
  rank-asymmetric multistart state that could hang a later collective
  bandwidth search without adding a separate synchronization protocol. Calls
  made before pool startup remain local, while an explicit
  `mpi.bcast.cmd(npseed(...))` continues to seed once per executing rank
  without a nested broadcast.

* Corrected active-MPI conditional-density derivative hat evaluation so its
  internal kernel-sum calls carry concrete argument values across the
  autodispatch boundary. This removes a rank-local lexical lookup failure
  while preserving the serial arithmetic, warnings, and chunking contract.

* Corrected active-MPI derivative kernel-weight output so every rank
  initializes non-owned cells before the canonical reduction. Output
  ownership is now invocation scoped and each derivative block retains its
  own tree-support intervals, also correcting compact-support adaptive-
  nearest-neighbour tree output and yielding exact local, one-worker, and
  multi-worker `p.kw` results without additional observation-square storage.

* Corrected positive-degree raw local-polynomial bandwidth selection for
  conditional density and distribution estimators when predictor units make
  the raw polynomial coordinates numerically ill conditioned (for example,
  calendar years). The signed full-row influence engine now selects one
  globally equivalent, conditioned coordinate representation before row
  traversal and uses it across scalar, block, tree, and all-large owners. The
  selection is basis-family neutral, does not center at evaluation points, and
  retains only `O(n k + k^2)` workspace for basis width `k`. The scalar row
  owner also now reuses the hoisted LP basis for the whole objective rather
  than rebuilding and clearing it for every evaluation row.

* Conditional count-bootstrap evaluation with a beta predictor or response
  kernel now uses one family-neutral native ingress for continuous and mixed
  categorical data, with strict `np.categorical.compress = TRUE`/`FALSE`
  behavior. Beta X or Y rows retain the canonical scaled-row owner, while the
  other side uses the shared continuous/categorical scalar registry; the
  implementation streams `O(n + B)` scratch and supports an entirely
  categorical non-beta side. For generalized- and adaptive-nearest-neighbour
  bandwidths, compressed multiplicities are expanded before row evaluation so
  bandwidth radii belong to the realized resample rather than the unreplicated
  source rows.

* Corrected conditional plot/bootstrap dispatch so local-polynomial beta fits
  retain the declared LP engine on every exact resample instead of being
  intercepted by the local-constant count evaluator. MPI exact refits execute
  locally within each outer bootstrap worker, avoiding nested collectives.

* Conditional count-bootstrap levels with a beta predictor or response kernel
  now reuse the canonical scaled conditional rows and contract all count
  columns through BLAS. The evaluator streams one observation-length row per
  side and retains no observation-by-evaluation matrix. This preserves the
  established weighted conditional ratio up to floating-point accumulation
  order while substantially reducing repeated count-evaluation time.

* Conditional-distribution least-squares bandwidth selection with a beta
  predictor or response kernel now uses the shared conditional
  local-polynomial CVLS engine across scalar and positive degrees, raw and
  Bernstein bases, fixed, generalized-nearest-neighbour, and
  adaptive-nearest-neighbour bandwidths, and supported mixed
  continuous/categorical data. Predictor rows use the canonical signed
  delete-one influence solve, while response rows use the canonical CDF
  operator on the actual evaluation plane. A bounded response-row supertile
  retains `O(n B + B^2)` workspace and reuses each CDF tile across small
  groups of predictor tiles; allocation failure selects only the
  algebraically identical linear-memory row stream. The former private beta
  conditional-distribution objective and its caller-free cache helpers have
  been removed.

* Conditional-density least-squares bandwidth selection with a beta predictor
  or response kernel now uses the shared conditional local-polynomial CVLS
  engine across scalar and positive degrees, raw and Bernstein bases, fixed,
  generalized-nearest-neighbour, and adaptive-nearest-neighbour bandwidths,
  and mixed continuous/categorical data. Predictor rows use the canonical
  signed full-row solve and exact delete-one identity; response rows retain the
  established bounded-quadrature or analytic-convolution owners. A bounded
  analytic-response supertile reuses convolution rows without retaining an
  observation-square matrix. Automatic Powell, MADS plus Powell, and degree
  search now accept these supported conditional-density configurations;
  conditional-distribution bandwidth activation remains separately gated.

* Conditional density and distribution gradients with a beta response kernel
  now retain the same signed-infinite endpoint limits and explicit missing
  gradient standard errors as the canonical beta regression engine. Finite
  gradients still use the strict finite restoration path, and indeterminate
  extended-real scaling remains an error; plot geometry continues to omit
  non-finite coordinates without mutating returned estimator values.

* Corrected internally computed generalized-nearest-neighbour convolution
  kernel sums so every continuous coordinate owns the required
  training-point companion bandwidths. The public route now uses the same
  unequal-bandwidth convolution contract as supplied-bandwidth consumers
  instead of dereferencing an absent companion matrix; fixed and adaptive
  bandwidth routes are unchanged.

* Corrected conditional-density least-squares cross-validation so its
  integrated-square and fitted-value terms use the same signed delete-one X
  smoother in every analytic and categorical-profile route. Adaptive-nearest-
  neighbour conditional rows now apply their observation-specific X
  bandwidth divisors and pair response convolutions with the corresponding
  row bandwidth. A scalar response-row shortcut is now restricted to ordinary
  density evaluation rather than being reused for convolution. These changes
  intentionally correct affected CVLS objectives and adaptive-nearest-
  neighbour CVML and conditional-distribution objectives; bounded-quadrature
  routes and adjacent fixed/generalized-nearest-neighbour objectives are
  unchanged.

* Beta density bandwidth selection now supports mixed continuous, unordered,
  and ordered data for CVML and bounded target-quadrature CVLS, across fixed,
  generalized-nearest-neighbour, and adaptive-nearest-neighbour bandwidths and
  beta orders 2, 4, 6, and 8. The continuous beta row and categorical kernels
  use the canonical density objective owners, including strict
  `np.categorical.compress = TRUE`/`FALSE` behavior; no beta-specific
  objective or categorical fallback is used. Other automatic beta bandwidth
  families retain their existing continuous-only gates.

* Corrected native automatic-degree search for conditional density and
  distribution bandwidths so every callback selects the canonical scalar
  degree-zero or general positive-degree local-polynomial engine from that
  callback's degree. Returned bandwidth objects now reproduce their stored
  objectives at the selected degree. Least-squares quantile degree search now
  updates public and engine metadata atomically, including singleton and
  Powell-refinement routes, instead of producing an inadmissible first
  objective from mixed legacy state. Automatic conditional-density callbacks
  also use a rank-invariant transition-consensus schedule.

* Conditional estimators and their gradient, hat, plotting, quantile, mode,
  and progress helpers now read canonical local-polynomial engine metadata by
  exact field name and validate the fields jointly. Malformed or mutated
  bandwidth objects fail with an error identifying the invalid field instead
  of reconstructing another engine from public display metadata. Derived
  least-squares quantile objects now retain the canonical engine state of
  their regression bandwidths.

* `np.categorical.compress` now requires a single non-missing logical value
  whenever categorical compression is applicable. Invalid explicit values
  fail with a precise option error instead of silently disabling the selected
  compression route; valid `TRUE`/`FALSE` behavior and routes where the option
  is inapplicable are unchanged.

* Fixed-bandwidth categorical conditional-density CVLS now bounds its three
  predictor, response, and response-convolution profile matrices within one
  checked workspace: 64 MiB for rank-one execution and 256 MiB for active
  MPI. Above the ceiling, stable resident rows and dedicated scratch rows are
  filled through the shared categorical-profile kernel engine. A 32-profile
  convolution supertile preserves every scalar accumulation order while
  reusing each response-convolution row across the group, reducing
  representative objective time without changing objective bytes. Rank-local
  readiness is agreed before collective traversal, and internal failure after
  route commitment is terminal rather than a silent fallback.

* Fixed-bandwidth categorical conditional-distribution CVLS now bounds its
  aggregate predictor- and response-profile kernel workspace at 64 MiB for
  serial/rank-one execution and 256 MiB under active MPI. Problems that fit
  retain the established dense arithmetic exactly; larger profile topologies
  use stable reuse-priority resident rows plus one scratch row, filled through
  the shared categorical-profile kernel engine. Objective accumulation order
  and values are unchanged, and route readiness and committed failure are
  rank-symmetric.

* Fixed-bandwidth ordered-categorical unconditional-distribution CVLS now
  consumes compressed kernel rows through a checked, caller-owned tile capped
  at 64 MiB for serial/rank-one execution and 256 MiB under active MPI instead
  of retaining the full evaluation-by-training profile matrix. The traversal
  and loss accumulation orders are unchanged, and objective values remain
  byte-identical; the general route remains canonical when a single row cannot
  fit the bounded workspace. Once the bounded route is selected, internal
  failure now fails rank-symmetrically rather than silently re-entering the
  general route.

* Autodispatch result sanitation now preserves explicit `NULL` call arguments
  without deleting their list slots while replacing temporary symbols. This
  corrects a shared post-computation failure that affected calls such as
  `npudistbw(..., gdat = NULL, do.full.integral = TRUE)` and could shift or
  overrun later arguments while reconstructing the returned call. Estimator
  arithmetic and MPI work ownership are unchanged.

* Added a dormant, rank-local categorical-profile tile engine for subsequent
  migration of dense profile consumers. Its caller-owned output is bounded by
  a checked 64 MiB ceiling, immutable training profiles can be validated once
  per traversal, and no production estimator dispatch changes in this
  checkpoint. The same work corrects unordered derivative dispatch so score
  operators reach their canonical kernels instead of the historical normal-
  kernel fallback; ordinary objective, fit, and prediction results and MPI
  ownership are unchanged.

* Rank-local Powell-side objective caches now use checked capacity, key-width,
  load-factor, and rehash arithmetic under an independent 64 MiB peak ceiling
  per native table. Growth is transactional on each rank, and the existing
  all-rank cache-hit agreement prevents a candidate cached on only some ranks
  from bypassing evaluation. Allocation-failure diagnostics are aggregated
  once when the search returns, avoiding new hot-loop synchronization. Cache
  keys, objective values, optimizer behavior, and the strict
  `np.objective.cache` option contract are unchanged.

* Smooth-coefficient local-polynomial cross-validation and fitting now solve
  common stable zero-ridge row systems through one registered rank-local
  native entry, reusing a bounded `O(p^2 + p)` LAPACK workspace instead of
  making one R-to-LAPACK transition per owned evaluation row. Successful
  wider batches also project their hoisted basis rows in ascending basis order
  through one register-local native entry instead of slicing, reshaping, and
  dispatching every row in R. Width-one and width-two explicit solvers are
  unchanged; if any wider row in a local batch is non-finite, singular, or
  ill-conditioned, that complete batch returns to the established R ridge
  loop. Estimator formulas, ridge policy, rank ownership, and collectives are
  unchanged, although LAPACK workspace alignment and projection reassociation
  can change last floating-point bits.

* Completed the migration to one canonical, uncentered local-polynomial
  compute engine across cross-validation, fitting, prediction/evaluation,
  gradients, standard errors, and hat/apply helpers. Raw and Bernstein bases
  now use the same basis-neutral solve machinery; the obsolete centered and
  Numerical Recipes matrix engines, including `linalg.c`, have been removed
  in favor of bounded-workspace BLAS/LAPACK implementations.

* Local-constant smoothing and explicit all-zero-degree local-polynomial
  smoothing now share the canonical implicit width-one LP engine throughout
  conditional objectives, fitting/evaluation, derivatives, hats, and helper
  routes. The scalar specialization constructs no unit basis and performs no
  BLAS/LAPACK solve; categorical-only routes continue to honor
  `np.categorical.compress = TRUE` and `FALSE`. Public regression-type
  metadata and estimator definitions are unchanged. This also corrects
  adaptive-nearest-neighbor conditional-density CVML for explicit degree zero:
  its historical general-LP row omitted observation-specific X-bandwidth
  divisors and could differ from the equivalent local-constant objective.

* Corrected conditional-density CVML workspace ownership for categorical
  responses when the scalar streaming objective is selected with multiple
  continuous predictors. Affected bandwidth searches, including
  `npconmode()` workflows, could previously abort the R process. Native MADS
  callbacks for the same scalar route also no longer allocate or read a
  nonexistent local-polynomial degree vector.

* Generic local-polynomial CV objectives now pack their immutable
  response-plus-basis operand once per rank and objective evaluation for
  eligible full drop-one rows, instead of repacking it for every owned
  observation. Adaptive-nearest-neighbour, tree, reduced-row, and specialized
  resident-row routes are unchanged, as are objective values.

* Conditional-density and conditional-distribution local-polynomial
  cross-validation now retain signed higher-order predictor-kernel weights in
  every delete-one X row on every rank. The obsolete QR route silently
  discarded negative weights; it has been removed in favor of the canonical
  signed full-row solve and exact diagonal deletion, with a sign-preserving
  denominator. This intentionally corrects affected higher-order CVLS and
  CVML objectives; ordinary second-order results and MPI ownership are
  unchanged.

* Conditional local-polynomial full-row deletion now preserves every finite,
  nonzero signed `1 - H_ii` denominator exactly rather than replacing
  sub-machine-epsilon values by a fixed floor. Exactly zero or non-finite
  denominators fail through the existing clean objective path because the
  corresponding deleted system has no valid finite row. Ordinary
  non-singular objectives and MPI ownership are unchanged.

* Eligible fixed and generalized-nearest-neighbour conditional
  local-polynomial full-row blocks now reuse the bounded weighted-design
  BLAS assembly already used by conditional-density CVLS. This accelerates
  conditional-density likelihood cross-validation and conditional-
  distribution least-squares cross-validation without changing the scalar
  width-one, adaptive-nearest-neighbour, non-Accelerate, or allocation-
  fallback paths. The formulas and MPI ownership are unchanged, but BLAS
  reassociation can change last floating-point bits.

* Adaptive-nearest-neighbour conditional local-polynomial objectives now
  reuse the same bounded weighted-design BLAS assembly on each rank when the
  basis has at least four terms. Widths one through three retain their
  established scalar transcript, as do small-sample, non-Accelerate, and
  allocation-fallback routes. This accelerates conditional-density likelihood
  and least-squares cross-validation and conditional-distribution
  least-squares cross-validation; formulas, rank ownership, and collectives
  are unchanged, but BLAS reassociation can change last floating-point bits.

* Adaptive-nearest-neighbour conditional-density likelihood
  cross-validation now reuses its rank-local X- and Y-row contexts across
  observations, matching the existing adaptive least-squares architecture.
  This removes repeated bandwidth and workspace setup without changing row
  ownership, reductions, or estimator formulas.

* Eligible adaptive-nearest-neighbour conditional-density objectives now
  cache ordinary-Gaussian response-bandwidth reciprocals in each rank-local
  row context, replacing repeated vector divisions with multiplications.
  Non-Gaussian, higher-order, bounded, mixed-response, tree, non-Accelerate,
  and allocation-fallback routes retain the established implementation.
  Formulas, rank ownership, and collectives are unchanged, but reciprocal
  multiplication can change last floating-point bits.

* On Apple ARM64, eligible fixed and generalized-nearest-neighbour products of
  two or more ordinary Gaussian continuous kernels now use one vector
  exponential of the summed squared standardized distances on each rank.
  Non-Apple, one-dimensional, tree, bounded, higher-order,
  adaptive-nearest-neighbour, convolution, score, permutation, and
  generalized-NN large-bandwidth shortcut routes retain the established
  implementation. The formulas are mathematically identical but can differ
  in their last floating-point bits.

* Fixed-bandwidth fourth- and sixth-order Gaussian convolution rows now hoist
  bandwidth- and evaluation-point-invariant polynomial terms outside their
  observation loop. The isolated helper preserves the established formulas,
  dimension-product order, and objective bytes; ordinary Gaussian,
  nearest-neighbour, bounded, score, permutation, and non-convolution routes
  retain their existing direct implementations.

* Corrected canonical local-polynomial CVAIC diagonal restoration and
  assembly for degrees greater than one, generalized-nearest-neighbor
  degree-one fitting and standard errors, and generalized-nearest-neighbor
  conditional-density CVLS convolution-bandwidth ownership. These are
  intentional numerical corrections and can change affected 0.70-6 results
  relative to 0.70-5.

* Local-polynomial solve retries are now bounded across objective, fit,
  prediction/evaluation, and hat-matrix paths. Non-finite or unrecoverable
  systems fail collectively instead of retrying indefinitely, and the fit
  route releases its native workspace before reporting failure; ordinary
  successful systems retain their established solve, ridge sequence, and MPI
  ownership.

* Aligned the dedicated `plot()` methods for `npregiv` and
  `npregivderiv` with the package-wide regression plotting vocabulary. They
  now use `gradients`, `data_overlay`, and `data_rug`; `npregiv` defaults to
  its structural level with the training response overlaid, while
  `npregivderiv` continues to default to its derivative. Automatic limits
  include every active curve, response overlay, and training-support rug,
  without replacing available fit-time evaluation curves. The experimental
  plot-only controls `plot.data`, `deriv`, and `phi` have been removed and now
  fail with migration guidance. Plotting remains master-side; estimator
  objects, MPI execution, and numerical results are unchanged.

* Conditional density and distribution level plots now support bootstrap
  intervals when the continuous X and Y sides use different kernel families
  or orders, including beta on either side. Beta-containing bootstrap levels
  use the same signed log-domain numerator/denominator accumulation as the
  fitted estimator; this also corrects same-family higher-order beta
  bootstrap values when a finite negative or extremely small explanatory
  kernel sum was previously replaced by a positive machine-epsilon guard.
  Fixed, frozen generalized/adaptive nearest-neighbour, and exact
  nearest-neighbour refit semantics are covered. Bootstrap replicates remain
  partitioned across MPI workers, and each worker evaluates its assigned
  counts locally without nested MPI dispatch. Matching legacy-kernel paths,
  estimator fits, gradients, bandwidth selection, and public defaults are
  unchanged.

* Aligned `npregiv()` and `npregivderiv()` summary output with `npreg()`.
  Fixed smoothing now uses the canonical `Kernel Regression Estimator:`
  description (`Local-Constant`, `Local-Linear`, or the local-polynomial
  degree and basis details) instead of separately printing the legacy `p`,
  regression type, and degree fields. This is a reporting-only change;
  fitted objects, smoothing choices, and numerical results are unchanged.

* MPI result simplification now applies one explicit scalar-logical contract
  at master-side ingress. Vector, empty, and missing `simplify` values are
  rejected before remote execution, parallel apply/replication, or simulation
  work is dispatched, preventing rank divergence and an unusable worker pool.
  Scalar logical values and numeric `0`/`1` retain their established results
  and MPI transport paths.

* `npksum()` now validates its public logical controls and `kernel.pow` before
  MPI or native dispatch. This prevents malformed scalar controls from being
  masked by short-circuit evaluation and prevents empty, missing, non-finite,
  or vector-valued kernel powers from reaching C; an empty kernel power could
  previously terminate the R process. Documented scalar inputs and valid
  integer powers retain their established behavior.

* Repaired `npksum()` score and OCG output for mixed data. `compute.score=TRUE`
  now allocates its native result buffer, and packed score, OCG, and continuous
  permutation blocks (including derivative kernel weights) are returned in
  original data-column order. Categorical base-kernel factors are now retained
  in mixed-data continuous permutation blocks, and score-only ordered kernels
  no longer dereference OCG-only state. MPI ranks now also assemble their
  disjoint derivative kernel-weight slices before returning `p.kw`, including
  uneven fixed/generalized-NN partitions and adaptive-NN training-row
  partitions. Uneven fixed/generalized-NN partitions now also use exact
  per-rank receive counts, preventing padding from one packed `p.ksum` block
  from overwriting the first element of the next. The defects could previously
  cause a process-level crash, fail during R reconstruction, or return an
  incomplete mixed-data derivative or derivative-weight array.

* Corrected Landweber--Fridman state coherence in `npregiv()`: state `N`
  now means exactly `N` completed updates, and `phi.mat[, N]`,
  `norm.stop[N]`, `norm.index`, the returned curve, derivatives, and optional
  weights refer to that same state. This can change the selected iteration and
  returned estimate. Bandwidth replay is now boundary-safe at `norm.index = 1`
  and retains recomputed stopping diagnostics. Documented numeric
  `starting.values` now initialize a complete level-and-derivative state.

* Corrected multivariate `npregiv()` derivatives to return one named pure
  coordinate partial per continuous structural coordinate. Local-constant
  Tikhonov derivative weights now use the derivative operator and bandwidth
  for the requested coordinate instead of combining coordinates or recycling
  bandwidth divisors. Univariate Tikhonov results are unchanged, and serial,
  one-worker, and multi-worker routes share the corrected contract.

* Separated `npregivderiv()` training and evaluation state. A fit-time
  `zeval`/formula `newdata` grid now affects only `phi.prime.eval` and its path;
  training operators, residuals, centering, stopping, selected state, fitted
  values, and training gradients remain unchanged. Arbitrary positive grid
  sizes are supported, while a different `weval` is rejected because the
  inverse problem is defined on the training instruments. The default Issue
  57 local-linear training route is numerically unchanged.

* Hardened the IV public surface in parity with `np`: ordinary IV regression
  stages now own `bandwidth.divide`, `ukertype`, and `okertype` without
  duplicate-argument failures and use kernel-appropriate categorical search
  bounds; derivative seeds, scalar controls, summaries, accessors, and plots
  now follow the documented contracts. MPI autodispatch now fingerprints all
  tagged remote arguments before reusing a worker reference, so modifying a
  returned bandwidth or IV replay object locally cannot silently reuse stale
  worker state. Post-fit `predict()` remains deliberately unavailable.

* The omitted regression-smoothing choice in `npregivderiv()` is now local
  linear (`regtype = "ll"`, degree one), matching the longstanding `p = 1`
  default of `npregiv()`. Explicit `regtype = "lc"` reproduces the former
  omitted-default computation, while `regtype = "lp", degree = ...` remains
  available. Local-polynomial order applies only to continuous predictors;
  categorical-only internal stages use the equivalent local-constant route.
  This intentional default change can alter derivative paths, stopping states,
  and fitted structural functions for calls that previously omitted
  `regtype`.

* `npregivderiv()` now forwards user-supplied unordered and ordered
  categorical regression kernels to its internal regression stages without
  colliding with the private Equation (14) adjoint. The adjoint continues to
  own its required Li--Racine categorical kernels and ordinary-CDF
  normalization; previously successful calls and estimator defaults are
  unchanged.

* `npregivderiv()` progress now identifies the mathematical object being
  computed, including `E[y|w]`, initialization derivatives, conditional
  residuals, and the derivative adjoint `T*{E[y-phi(z)|w]}`. State-zero work
  is distinguished from work at completed iteration `N`; estimator results,
  stopping behavior, smoothing choices, and MPI lifecycle are unchanged.

* Modernized the public `npregiv()` and `npregivderiv()` interfaces in parity
  with serial `np` while preserving their established numerical defaults.
  Both now support explicit IV formulas (`y ~ z | w`, with optional `| x`
  where the estimator supports it), `data`, `subset`, `na.action`, fit-time
  `newdata`, fixed `regtype = "lc"`/`"ll"`/`"lp"` and scalar `degree`
  controls, structured summaries, and `fitted()`, `gradients()`, and
  training-row `residuals()`. Objects retain compact bandwidth, smoothing,
  and stage metadata, and the formula calls remain rank-symmetric under MPI
  autodispatch. Legacy `p` remains supported by `npregiv()` with unchanged
  default behavior. Unsupported post-fit prediction, derivative `x`, and
  automatic-degree NOMAD routes are documented and fail explicitly instead
  of being silently accepted or approximated.

* Beta associated kernels now interpret \code{ckerbound="range"} using
  outer half-spacing bounds based on the two smallest and two largest distinct
  training values. This keeps every observation strictly inside the fitted
  support, removes empirical-extremum jumps and infinite raw-extremum
  derivatives, and preserves tied-extremum multiplicity. Explicit fixed beta
  bounds remain literal, while every non-beta range route continues to use
  exact sample minima and maxima.

* Beta regression hat matrices, leverages, and matrix-free applications now
  use the same canonical common-scaled row and local-polynomial solve owners as
  fitted values. Orders 2/4/6/8, fixed and nearest-neighbour bandwidths,
  raw/Bernstein bases, derivatives, and mixed categorical predictors are
  supported without reconstructing underflowed absolute weights. Direct
  application retains only a bounded influence-row block rather than an
  evaluation-by-training matrix; a full matrix is allocated only when the
  caller explicitly requests one.

* Beta range resolution now validates distinct-extrema metadata coordinate by
  coordinate. This fixes an unreleased 0.70-6 development defect in which two
  or more continuous beta coordinates failed during support resolution because
  vector and scalar logical operators were mixed.

* Automatic univariate second-order beta range searches now certify fixed-
  bandwidth density CVLS and distribution CDF solutions against the
  double-precision uniform-limit objective. A material improvement triggers
  one deterministic Powell refinement from the resolved support width, and
  the best material candidate is retained. The user's ordinary search and
  \code{nmulti} setting are unchanged, and cold starts and restarts from
  bandwidth objects share the same contract; summaries report the additional
  certification evaluations. Explicit fixed bounds, CVML, higher beta orders,
  nearest-neighbour modes, multivariate searches, and non-Powell solvers are
  unchanged.

* Corrected recursive fitted-value centering in `npregivderiv()` so both
  empirical terms of the Equation (14) adjoint use the same fitted
  conditional-residual vector. The recursive path had instead centered its
  second term on the raw residual even though its first term used the fitted
  residual. Serial and MPI autodispatch routes share the corrected contract.
  The `npregiv()` and `npregivderiv()` examples now both use `n = 500`, for
  which the documented seed produces stable estimates.

* Corrected Landweber-Fridman state indexing in `npregivderiv()`. Iteration
  `N` now consistently denotes `N` completed derivative updates: column `N`
  of `phi.prime.mat` and `phi.mat`, `norm.stop[N]`, `num.iterations`, and the
  returned derivative/function now identify the same state. The initialization
  remains state zero, and evaluated stopping-rule overshoots remain available
  in the iteration matrices while `num.iterations` identifies the selected
  state. This can change the selected iteration and estimates because the
  stopping criterion is now evaluated as `N` times the residual norm at state
  `N`, as described by Florens, Centorrino, and Racine. Serial and MPI
  autodispatch routes share this contract.

* Corrected the empirical adjoint in `npregivderiv()` to use the ordinary
  kernel CDF required by Equation (14) of Florens, Centorrino, and Racine.
  The integral kernel sum had retained an extra continuous-bandwidth factor,
  which could make the Landweber-Fridman stopping norm increase and bias the
  recovered conditional mean. The private adjoint now owns its required
  normalization, while public `npksum()` defaults, regression argument
  forwarding, and MPI payloads remain unchanged. The iteration guard also
  examines only the computed prefix of the preallocated stopping vector. This
  mirrors the serial `np` repair and resolves issue #57.

* Retired the unused experimental truncated-Gaussian continuous kernel and
  its `nptgauss()` configuration helper. The supported continuous kernels are
  Gaussian, Epanechnikov, and uniform; their public interfaces, native codes,
  numerical behavior, and MPI payloads are unchanged.

* Corrected heterogeneous generalized local-polynomial (GLP) construction to
  include every coordinate-capped term through the declared total degree. For
  example, `degree = c(2, 1)` now includes the `x1*x2` term. Raw and
  `bernstein.basis = TRUE` GLP fits now use exact representations of the same
  complete polynomial space; the latter uses a deterministic degree-graded,
  orthonormal shifted-Legendre representation. Raw heterogeneous GLP results
  can therefore differ from earlier versions when terms are restored;
  multivariate `bernstein.basis = TRUE` GLP and automatic degree-search results
  can also differ because the former fixed-degree Bernstein columns did not
  span the corrected GLP space once interactions entered. Univariate and
  maximum-degree-one GLP retain the same polynomial space, although their
  floating-point paths can differ. Dimension guards, higher-order
  derivatives, fitted values, standard errors, and the serial/MPI native
  evaluators use the same term definition. Local constant, local linear,
  univariate, additive, tensor, and non-local-polynomial semantics are
  unchanged.

* Hardened native NOMAD observer and interrupt handling. Explicit user
  interruption is now reported as an R `interrupt` condition only after
  native cleanup; active MPI work defers a master interrupt until the current
  rank-common computation boundary so workers remain reusable without adding
  MPI commands or payloads. Ordinary observer errors remain fail-open, native
  callables are resolved per solve rather than retained across package
  reloads, and the declared R 3.5 compatibility floor is preserved.

* Restored timely master-only progress updates during native NOMAD bandwidth
  and degree searches, including long compiled objective evaluations.
  Iteration, current degree, and accepted-best details now follow the
  package-wide progress interval without changing optimization results,
  evaluation accounting, MPI commands, or payloads. This requires `crs`
  0.15-46 or later.

* Restored partial mixed-degree local-polynomial gradient evaluation and
  plotting across regression, conditional density/distribution, least-squares
  quantile regression, conditional quantiles, and conditional modes. A
  continuous derivative requested above that coordinate's fitted degree is
  retained as `NA` and shown as an empty panel while other supported
  derivatives and categorical first-difference effects continue to plot.
  Requests with no available component now fail before bootstrap or MPI work,
  with the affected predictor orders and degrees and actionable refitting
  guidance. MPI commands and payloads are unchanged.

* Corrected common-scale uncertainty-band rendering for partially available
  gradient plots. Empty derivative panels now retain their list positions and
  finite-range calculations no longer emit spurious `min`/`max` warnings;
  supported estimates, interval bands, and MPI payloads are unchanged.

# npRmpi 0.70-5

* Corrected the heteroskedasticity-robust Ichimura index-coefficient
  covariance returned by `vcov(npindex(..., gradients = TRUE))` when the model
  has more than one free index coefficient. Link-gradient and residual weights
  are now applied observation by observation. Fits with one free coefficient
  remain exactly unchanged, as do beta, bandwidth, objective, fitted-value,
  residual, and gradient results. Serial and MPI results remain exactly
  equivalent, with no payload change.

* Repaired finite-support bounded continuous-kernel normalization for all
  supported Gaussian, Epanechnikov, uniform, and truncated-Gaussian kernels.
  The centered normalization now avoids tail cancellation, uses the platform
  C99 error function where appropriate, and preserves the uniform-density
  limit as the bandwidth tends to infinity.

* Replaced the bounded distribution operator with the mathematically correct
  observation-centered truncated-kernel CDF. Analytic centered interval
  primitives and hoisted invariant work preserve exact support endpoints,
  numerical accuracy, the finite-support large-bandwidth uniform limit, and
  serial/MPI numerical equivalence.

* Corrected `transform.bounds = TRUE` initialization for regression initial
  starts and multistarts and for later conditional-distribution multistarts.
  External bandwidth starts are now inverse-mapped before transformed Powell
  search; the public default remains `FALSE` and MPI payloads are unchanged.

* Clarified optimizer summaries so objective-cache lookups, NOMAD point
  lookups, family-native R-`optim` refinement, and total function-evaluation
  accounting are labelled consistently without implying that their
  denominators are interchangeable.

* Clarified bandwidth-object replay documentation across estimator families.
  Stored bandwidth objects retain search and fit metadata but do not silently
  materialize or expand training data; callers must provide data again where
  the documented replay route requires it.

* Hardened MPI lifecycle and protocol handling. Bootstrap fanout timeout or
  failure now reports a hard failure and instructs the user to restart R before
  further MPI-backed computation; attach close resets pool state, profile
  communicator metadata remains distinct, and internal protocol tags are
  allocated portably.

* Standardized the public `nomad.opts` contract across supported bandwidth
  selectors and hardened native NOMAD callback cleanup and allocation guards
  without changing MPI message payloads.

* Refreshed MPICH build and runtime guidance for current MacPorts layouts and
  documented the supported session, attach, and profile/manual-broadcast
  launch modes. Open MPI remains outside the validated 0.70-5 backend matrix.

* Reworked `npindexbw()` / `npindex()` internals after 0.70-4. Ichimura and
  Klein-Spady single-index objectives now reuse the established `npreg`
  leave-one-out backend where applicable, preserving the public estimator
  contract while materially improving high-dimensional and local-polynomial
  objective evaluation. MPI-aware objective services now keep the root-owned
  optimizer and worker evaluation protocols explicit.

* Repaired MPI single-index NOMAD/Powell handoff behavior for
  `npindex(..., nomad = TRUE)` and related bandwidth searches. The Powell
  refinement stage now uses the same one-active-protocol discipline as the
  validated MPI search routes, avoiding worker/master desynchronization during
  the NOMAD-to-Powell transition.

* Distributed post-CV `npindex()` fitting, evaluation, gradients, and related
  components across owned evaluation rows in active MPI sessions while keeping
  the full training data available where kernel smoothing requires it. Focused
  large-sample checks showed the expected scaling improvement and preserved
  serial-equivalent numerical results.

* Added user documentation for `npindex()` optimizer choice and beta
  interpretation. The documentation now gives practical guidance on when
  derivative-free Nelder-Mead remains a reasonable low-dimensional default and
  when BFGS is useful for higher-dimensional index searches, and it explains
  the relative interpretation of normalized single-index beta coefficients.

* Hardened single-index fit, evaluation, plotting, service, and summary
  behavior. Bounded continuous-kernel options are now carried consistently
  through objective, fit, evaluation, variance, and bootstrap routes;
  Klein-Spady confusion-matrix output is guarded against out-of-range fitted
  values; and large-bandwidth shortcuts are restricted to kernel/order
  combinations whose constant-weight approximation was validated.

* Repaired single-index plot bootstrap memory usage. Wild-bootstrap
  self-maps now avoid public-facing `O(n^2)` kernel-weight allocations for the
  ordinary large-sample plot route, while preserving public fitted/plotting
  contracts and MPI row ownership.

* Repaired generalized-nearest-neighbor local-polynomial derivative ownership
  across `npreghat()` and related public routes. Mixed-degree local-polynomial
  fits now route available derivative components through the correct owner,
  preserve `apply == H %*% y` contracts, and report unavailable derivative
  components consistently rather than silently applying the wrong operator.

* Repaired local-constant derivative ownership in `npreghat()` so degree-zero
  local-polynomial derivative requests use the analytic local-constant
  derivative contract only when that contract is mathematically available.
  Scalar, matrix, multi-column, serial-equivalent, and active-pool error routes
  were validated separately.

* Improved conditional mixed-degree local-polynomial gradients for
  `npcdens()` / `npcdist()` and public conditional-gradient accessors. The
  partial-availability contract now matches the regression-family policy:
  available components are returned, unavailable components are represented as
  `NA`, and incoherent metadata fails clearly.

* Repaired MPI conditional-density native derivative workspace handling so
  no-gradient and gradient calls can be mixed safely in one session without
  stale native derivative workspace state affecting later calls.

* Repaired exact nearest-neighbor bootstrap parity in MPI conditional-density
  and conditional-distribution plot routes. Adaptive and generalized
  nearest-neighbor exact bootstrap paths now rebuild bandwidth state on the
  expanded resample consistently with the serial package.

* Repaired formula/data reentry contracts across density, distribution,
  conditional density, conditional distribution, single-index, smooth
  coefficient, partially linear, quantile, conditional-mode, copula, and
  significance-test routes. Explicit estimator `data=` now overrides stored
  bandwidth-object data where that public call shape is supported, formula
  `newdata` is validated against fitted RHS variables, and direct formula
  calls with numeric smoothing parameters no longer misroute formula objects
  as native data before autodispatch.

* Hardened MPI autodispatch argument transport and bandwidth fingerprints.
  Named arguments materialized through S3 dispatch are preserved before worker
  fanout, and bandwidth-object reuse checks are stronger against stale remote
  object reuse.

* Repaired unconditional density and distribution edge contracts. Ordered
  kernel-code selection is now consistent between conditional bandwidth
  selection and fitting where normalization is required, all-NA input is
  rejected before native calls, training and evaluation omission metadata are
  retained separately where needed, categorical zero-bandwidth standard-error
  handling is consistent across categorical configurations, and unsupported
  bandwidth-selection method codes fail clearly.

* Clarified normal-reference bandwidth documentation for density and
  distribution routes. The rule-of-thumb formulas are documented as fast
  exploratory Silverman-style heuristics, not production substitutes for
  cross-validation or likelihood-based selection.

* Hardened `npcdistbw()` normal-reference method handling so the R method code
  and native C method code stay aligned, avoiding accidental fall-through to
  an unintended bandwidth-selection branch.

* Repaired and hardened public estimator contracts found during adversarial
  audits of `npscoef`, `npplreg`, `npqreg`, `nplsqreg`, `npconmode`,
  `npcopula`, `npreg`, `npcdens`, `npcdist`, `npudens`, and `npudist`.
  Repairs include `npscoef` iterated backfit behavior and MPI service-task
  handling, partially-linear fit reentry, quantile inversion/clamping
  contracts, conditional-mode probability/tie handling and orchestrator
  localization, copula sample reentry, regression tree/large-h predicate
  alignment, and density/distribution formula/native argument consistency.

* Plot-bootstrap memory hardening now covers single-index, conditional-mode,
  and partially-linear plot routes. The default plot evaluation grids remain
  linear in the training sample size, and explicit `neval == ntrain` style
  requests remain user-controlled.

* Native shadow-object and cache lifecycles were hardened. Conditional-density
  shadow pointers, native objective-cache state, and regression large-h /
  large-lambda caches are cleared at the appropriate top-level lifecycle
  boundary so pointer-keyed helper state cannot leak across independent calls,
  datasets, or MPI sessions.

* `npksum()` numeric-bandwidth dispatch now constructs data-aware bandwidth
  objects for mixed continuous/categorical data instead of falling back to a
  default-typed object. The serial and MPI packages now share the same public
  behavior for numeric smoothing parameters.

* The local-polynomial regression CVKS low-support objective path now uses a
  named bandwidth-method code and aligns the R and C routing contracts,
  reducing the risk of drift between search metadata and native objective
  handling.

* Expanded focused tests, demos, benchmarks, and release-protocol sentinels
  across the public exported surface. The release protocol now requires an
  explicit public-exported-surface inventory, estimator-family sentinels,
  documentation/demo/benchmark smoke coverage, MPI scaling canaries, and
  installed package proof before release-ready claims.

# npRmpi 0.70-4

* Hardened proactive C cleanup paths by aligning LAPACK helper out-of-memory
  cleanup with the serial package, clearing extended nearest-neighbor alias
  state in the shared estimator cleanup helper, routing selected
  bandwidth-constructor allocation and unsupported-method failures through
  existing cleanup labels, and guarding large `np_kernelsum()` allocation-size
  products before calling the legacy `alloc_vecd(int)` allocator.

* Added `nomad = "auto"` for local-polynomial degree searches, preserving the
  serial package policy while keeping MPI autodispatch and NOMAD shadow
  protocols explicit. Small one-dimensional degree lattices use exhaustive
  Powell-style search; larger or explicitly requested NOMAD surfaces continue
  to use NOMAD.

* Repaired `npRmpi` fit-wrapper NOMAD guards after normalized NOMAD control
  handling so `nomad = TRUE` and `nomad = "auto"` routes stay on the intended
  local/native path rather than being accidentally autodispatched.

* Bandwidth and fit summaries now report cumulative search diagnostics more
  clearly. NOMAD cache output distinguishes repeated point lookups avoided by
  NOMAD from true objective computations; Powell summaries expose repeated
  objective lookups avoided by the package-side cache; hybrid NOMAD+Powell and
  exhaustive Powell timing labels are reported consistently.

* Base-graphics plot scaling and legends have been hardened for multi-panel
  displays. Plots now honor active `mfrow`/`cex` behavior more consistently,
  use role-appropriate legend sizes, and draw factor legends with point glyphs
  matching the plotted estimates.

* Fixed-bandwidth bias-corrected bootstrap plot intervals have been reworked
  across the supported serial-equivalent and MPI-aware plot families.
  Bias-corrected centers and intervals now share a common centering contract
  for regression, unconditional density and distribution, conditional density
  and distribution, single-index, partially linear, and smooth-coefficient
  plots where supported.

* Pair/block/geometric bootstrap intervals in regression-style routes now use
  smooth-bootstrap bias correction when `center = "bias-corrected"` is
  requested; wild-bootstrap regression intervals retain the standard
  wild-bootstrap correction. Density and distribution routes use
  perturbation-based smooth-bootstrap bias correction. Gaussian, uniform, and
  second-order Epanechnikov perturbation kernels are supported; higher-order
  signed perturbation kernels fail closed.

* Bias-corrected plot support now covers mixed-data unconditional density and
  distribution, mixed-data conditional density and distribution, and
  conditional gradient displays. Non-fixed/adaptive/generalized bandwidth
  bias correction and empirically unsupported `npqreg` bias correction remain
  fail-closed with explicit messages.

* Bias-corrected plot payloads returned by `output = "data"` and
  `output = "both"` now consistently expose fitted values, bias-corrected
  values, gradients, gradient bias corrections, and interval payloads with the
  same centering contract used for rendering.

* Derivative-order validation has been tightened across plot/gradient routes
  so unsupported derivative requests fail early or are represented
  consistently rather than silently plotting a lower-order derivative.

* Single-index formula dispatch with explicit bandwidths has been repaired in
  MPI-aware routes, and single-index bias-corrected plot centers, output
  payloads, and legends have been tightened.

* The native `loadNamespace("crs")` call construction used by CRS-backed
  native search routes is now protected, resolving the RCHK protection finding
  without changing the native NOMAD search contract.

# npRmpi 0.70-3

* Added MPI-aware `nplsqreg()`/`nplsqregbw()` support for location-scale
  quantile regression, including formula/data and bandwidth-object workflows,
  scalar/vector `tau`, prediction, residual extraction, summaries, and plot
  routes built on the shared quantile plotting engine.

* Supported MPI MADS/NOMAD-backed bandwidth-search routes now use the final
  native `crs` NOMAD C API rather than the retired legacy `snomadr()`
  fallback. The runtime dependency on `crs` is now declared in `Imports`,
  while `LinkingTo` remains for the native header.

* Native NOMAD routes now preserve progress best-record reporting, expose
  cache/evaluation diagnostics, honor explicit start and option controls, and
  reject unsupported or indeterminate cache-off settings before solver entry.
  Inadmissible GLP degree candidates are guarded before expensive evaluator
  work in serial-equivalent and MPI-dispatched routes.

* `npindexbw(..., method = "ichimura", regtype = c("ll", "lp"))` now reuses
  the established local-polynomial regression objective evaluator, and MPI
  autodispatch uses a rank-0-driven objective service for the fixed-degree and
  NOMAD degree-search Ichimura local-polynomial routes. Focused sentinel runs
  preserved payloads while restoring useful scaling for the formerly flat
  local-linear and local-polynomial single-index rows.

* MPI fanout coverage has been extended for computationally heavy bootstrap
  workloads in specification, dependence, distribution-equality, quantile, and
  symmetry tests, and plot-bootstrap RNG streams now restore the
  serial-equivalent final state after MPI fanout.

* The shipped `npplreg` attach-mode demo now explicitly finalizes the master
  rank after successful attach shutdown, hardening release-sentinel teardown
  without changing estimator or runtime defaults.

* MPI auto-dispatch for `nplsqreg()` now materializes named method-level
  `...` arguments before worker fanout, preserving user-supplied `scale` and
  option values that arrive through S3 `..n` placeholders.

* `options(np.tree = "auto")` is now the default tree mode. In auto mode,
  continuous kd-tree routes are enabled only for bounded-support continuous
  kernels (`"epanechnikov"` and `"uniform"`); `np.tree = TRUE` remains the
  explicit force-on override and `np.tree = FALSE` remains the force-off
  diagnostic path.

* Powell bandwidth searches now expose package-side repeated-candidate
  objective caching through `options(np.objective.cache = TRUE/FALSE)`. The
  cache remains enabled by default and is scoped to one bandwidth solve, so it
  can reuse exact candidates across Powell restarts without carrying state
  across datasets or later calls. Continuous-only generalized/adaptive
  nearest-neighbor routes also retain their integer nearest-neighbor objective
  cache under the same switch. The option is synchronized to MPI workers in
  autodispatch sessions; NOMAD solver caching and extended-NN distance reuse
  remain separate mechanisms.

* Continuous large-bandwidth shortcut evaluations can now be disabled with
  `options(np.largeh = FALSE)`, and discrete near-upper-bandwidth shortcut
  evaluations can now be disabled with `options(np.largelambda = FALSE)`.
  Both remain enabled by default and are synchronized to MPI workers in
  autodispatch sessions. These switches are intended for diagnostic timing and
  reproducibility studies that need to separate tree effects from
  large-bandwidth and large-lambda fast paths without changing the canonical
  dense/tree objective machinery.

* Local-polynomial regression cross-validation now uses a leaner hot
  symmetric weighted-sum loop. Fixed-bandwidth `npregbw(..., regtype = "lp",
  bwmethod = "cv.ls")` objective probes in active MPI sessions match serial
  `np` objective values to numerical precision while substantially reducing
  local-polynomial CV evaluation time.

* Shared weighted outer-product accumulation in `npksum()` now uses a guarded
  BLAS `dgemm` route when the operation is dense, non-permuted, and
  memory-bounded. Focused fixed-bandwidth probes preserve serial/MPI objective
  parity while substantially accelerating high-basis local-polynomial
  regression and smooth-coefficient objective rows; small and scalar routes
  remain on the established loop path.

* Unconditional density least-squares cross-validation now uses a leaner
  fixed-bandwidth Gaussian convolution loop. Fixed-bandwidth
  `npudensbw(..., bwmethod = "cv.ls")` objective probes preserve objective
  values exactly in the focused validation rows while materially reducing the
  convolution portion of the objective calculation. Conditional-density
  least-squares objective probes inherit the same fixed-bandwidth Gaussian
  convolution improvement.

* Non-Gaussian scalar-bandwidth convolution helpers now hoist the response
  bandwidth power outside the inner loop, improving fixed-bandwidth
  least-squares density cross-validation with compact-support kernels while
  preserving objective values exactly in focused probes.

* Continuous-kernel vector helpers now reuse the loop-invariant signed inverse
  bandwidth scale inside their inner loops. Focused density, conditional
  density, and regression objective probes preserve serial/MPI objective
  parity while reducing repeated scaling work in shared C hot paths.

* Conditional density and conditional distribution least-squares
  cross-validation now use a size-aware row-block policy for local-polynomial
  objective evaluation. The accepted route keeps the bounded-quadrature cap
  unchanged, bounds transient memory by sample size, and preserves objective
  values to numerical precision while materially reducing evaluator overhead
  for fixed-bandwidth CVLS probes in serial and MPI sessions.

* Local-polynomial conditional density maximum-likelihood cross-validation now
  uses the same bounded-memory block machinery for fixed and generalized
  nearest-neighbor bandwidths. Focused `npcdensbw(..., bwmethod = "cv.ml",
  regtype = "lp")` probes preserve objective values and selected bandwidths to
  numerical precision in serial and MPI sessions while reducing objective and
  full-search runtime.

* Large-sample categorical-only regression now uses the MPI-safe
  profile-compressed route under `options(np.categorical.compress = TRUE)`,
  which is enabled by default. This categorical route is independent of
  `options(np.tree)`. Repeated predictor profiles are compressed before
  bandwidth search, fitting, prediction/evaluation, standard errors,
  hat-helper use, and plot bootstrap helpers, preserving the established
  dense-route numerical contract while reducing repeated work.

* Categorical-only unconditional density routes now use the same
  profile-compression idea when `options(np.categorical.compress = TRUE)` is
  enabled. The fixed-bandwidth fit/evaluation route preserves dense-route
  fitted/evaluation values while avoiding repeated computation over identical
  categorical profiles, and the bandwidth-search route now uses the same
  compressed support representation for all-categorical data. As with other
  flat categorical search surfaces, selected smoothing parameters may drift by
  optimizer-path amounts while preserving the objective scale. Very fast
  compressed routes may remain overhead-floor limited, so MPI acceleration is
  most useful once the uncompressed work would be genuinely long-running.

* Categorical-only conditional density and conditional distribution bandwidth
  searches now honor `options(np.categorical.compress = TRUE)`. The promoted
  route preserves the objective value to numerical precision while allowing
  harmless optimizer-path drift in selected smoothing parameters, especially
  near upper-bound or large-bandwidth regions where the objective is flat.

* Ordered-only unconditional distribution bandwidth search and fit/evaluation
  routes also use profile compression when
  `options(np.categorical.compress = TRUE)` is enabled. The bandwidth-search
  route preserves the objective value to numerical precision while allowing
  harmless optimizer-path drift in selected smoothing parameters; fitted
  distribution values and standard errors are preserved while avoiding repeated
  computation over identical ordered profiles.

* Fixed-bandwidth local-constant `npscoef()` fits now use categorical-profile
  compression when all `Z` variables are categorical and
  `options(np.categorical.compress = TRUE)` is enabled. The route preserves
  fitted means, coefficient surfaces, asymptotic mean standard errors, and
  coefficient/gradient standard errors for training and evaluation fits while
  avoiding repeated work over duplicate `Z` profiles. The corresponding
  `npscoefhat(output = "apply")` path and count-based plot-bootstrap helper
  use the same profile compression without changing the explicit full-matrix
  `output = "matrix"` contract.

* Internal categorical-profile and large-bandwidth caches are now cleared at
  the relevant top-level density, distribution, conditional-density,
  conditional-distribution, and regression cleanup points. These caches are
  keyed by call-local row pointers, so clearing them per `.Call` prevents stale
  same-process state from leaking across unrelated data sets or MPI dispatch
  modes.

* Fixed `npcdens()` and `npcdist()` formula calls with explicit numeric
  smoothing parameters, such as `npcdist(y ~ x, data = dat, bws = c(.25,.25))`,
  so `npRmpi` preserves the established formula-to-bandwidth-object rewrite
  before MPI autodispatch.

* Hardened the `npudist()` formula route so formula calls are handled before
  MPI autodispatch.

* `npplreg()` now activates the already validated categorical regression
  compression path for its internal all-categorical `Z` regressions when
  `options(np.categorical.compress = TRUE)` is enabled, without requiring users
  to request continuous kd-tree acceleration through `options(np.tree)`.

* Formula variables whose names contain dots, such as `y.irr ~ x`, are no
  longer mistaken for the formula wildcard `.` in conditional density and
  conditional distribution bandwidth routes. The conditional-density bandwidth
  formula route also now expands the actual wildcard form `y ~ .` using the
  supplied `data` frame, matching the conditional-distribution route.

* Fixed MPI conditional-density and conditional-distribution NOMAD degree-search
  routes so Powell refinement and promoted wrappers such as `npconmode()` reach
  the intended bandwidth-object construction path rather than the pre-search
  autodispatch preflight used by non-degree-search routes.

# npRmpi 0.70-2

* `npqreg()` is now a fully fledged MPI-aware quantile-regression front
  end. It supports the formula/data workflow, internally computes
  `npcdistbw()` bandwidths when a bandwidth object is not supplied,
  accepts scalar or vector `tau`, reuses selected bandwidths for
  additional quantiles in `plot()`, and exposes the usual S3 surface:
  `fitted()`, `predict()`, `predict(..., se.fit=TRUE)`, `se()`,
  `gradients()`, `summary()`, `print()`, `quantile()`, and `plot()`.

* `npqreg()` prediction now honors the standard `newdata` workflow while
  preserving native `exdat` precedence for compatibility with existing
  `npRmpi` call surfaces. Formula-based prediction validates that new
  data contain the required right-hand-side variables.

* `npqreg()` plotting has been expanded for vector quantiles,
  level/gradient displays, ordered predictors, user-specified legends,
  and object-fed plotting of additional `tau` values without recomputing
  cross-validation. The fixed-bandwidth gradient path now uses the
  MPI-aware helper route.

* `npconmode()` is now a first-class conditional-mode estimator. It
  supports formula/data and bandwidth-object workflows, forwards
  bandwidth-selection options to `npcdensbw()`, propagates local
  polynomial and NOMAD metadata, and exposes `fitted()`, `predict()`,
  `summary()`, `print()`, `gradients()`, and `plot()` methods.

* `npconmode()` now supports optional class-probability matrices and
  level-specific probability gradients. For non-local-constant fits,
  probabilities are normalized to be non-negative and to sum to one
  across the discrete response support before modal classification.

* `npconmode()` now fails early for non-categorical responses and
  validates formula-based `newdata` against the original right-hand-side
  variables.

* `npconmode()` plotting now supports object-fed class-probability slices
  and two-dimensional probability surfaces, optional `rgl` rendering, and
  probability-level asymptotic intervals where defined. Surface bootstrap
  intervals for class probabilities remain intentionally deferred.

* `npcopula()` is now a first-class copula estimator. It supports
  formula/data and bandwidth-object workflows, automatic two-dimensional
  probability grids, explicit `u` evaluation grids, and ordinary
  extractable object components including `$bws`.

* `npcopula()` now provides `fitted()`, `predict()`, `predict(...,
  se.fit=TRUE)`, `se()`, `summary()`, `print()`, `as.data.frame()`, and
  richer `plot()` methods. Plotting supports base `persp`, `image`, and
  optional `rgl` rendering, with asymptotic and MPI-fanned bootstrap
  intervals for copula surfaces where defined.

* `npcopula()` explicit-grid evaluation now uses the direct estimator
  route, preserving numerical results while avoiding the severe runtime
  growth of the previous expanded-grid path when users request larger
  probability grids.

* The automatic local-polynomial NOMAD controls have been split into
  explicit restart toggles: `powell.remin` for Powell restarts and
  `nomad.remin` for the second NOMAD hot start. This preserves the
  Powell Numerical Recipes restart default while allowing NOMAD hot
  starts to be controlled separately.

* Deprecated legacy `remin` remains accepted by `npregbw()` and `npreg()`
  with a warning and is mapped to the modern `powell.remin`/`nomad.remin`
  controls where appropriate, preserving downstream compatibility while
  documenting the new spelling.

* Hat-operator helpers now support an additional constraint-oriented
  output route for objects needed by shape-constrained quadratic
  programming workflows, avoiding reimplementation of local-polynomial
  hat-matrix construction in user examples.

* Local-polynomial derivative support has been broadened across the
  conditional estimator family. `npreg()`, `npcdens()`, and `npcdist()`
  now honor `gradient.order` more consistently for fitted, evaluated,
  predicted, and plotted objects when the selected polynomial degree is
  high enough, including vector derivative orders over continuous
  predictors and tensor/additive/Bernstein local-polynomial bases. The
  MPI implementation dispatches the corresponding conditional hat-apply
  helper work across the active worker pool where applicable.

* Core and semiparametric S3 prediction paths have been hardened around
  `newdata`, native evaluation-argument precedence, formula RHS
  validation, and `se.fit` handling while preserving `npRmpi` route
  independence.

* Front-end/bandwidth argument hygiene has been tightened so
  estimator-only controls such as `proper` are not forwarded into
  bandwidth selectors that do not accept them.

* MPI lifecycle and plotting routes received additional hardening,
  including soft `npRmpi.quit()` behavior, local object-fed plot
  computation where required, and explicit fanout of applicable
  bootstrap workloads.

* Documentation has been refreshed for the promoted `npqreg()`,
  `npconmode()`, and `npcopula()` workflows, including the
  local-polynomial NOMAD route, probability/gradient outputs, plot
  controls, and examples that use the streamlined interfaces.

* The pre-release validation suite was expanded with focused hostile
  argument tests, S3 contract tests, installed/tarball proof scripts,
  route-aware MPI probes, and serial/MPI parity checks for the newly
  promoted estimator families.

# npRmpi 0.70-1

* The default multistart cap for bandwidth selection now follows
  `min(2, p)` across the mirrored estimator families, replacing the
  older `min(5, p)` cap. This includes automatic LP degree-search calls
  when `search.engine="nomad"` or `"nomad+powell"` and `nmulti` is not
  supplied explicitly.

* The univariate boundary density helper `npuniden.boundary()` now
  defaults to `nmulti=1`.

* The empirical studies supporting this mirror change are documented in
  `np-master/benchmarks/validation/`, with a summary note kept in this
  repository's `benchmarks/validation/` folder.

* LP-capable front ends now accept `nomad=TRUE` as a documented
  convenience preset for the recommended automatic NOMAD
  local-polynomial route, mirroring the serial package defaults and
  help-page guidance.
