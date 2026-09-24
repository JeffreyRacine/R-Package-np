# np 0.80-1

* Generalized-neighbor conditional-distribution cross-validation on external
  grids now uses the deleted donor sample for both predictor and response
  radii, retaining the empirical-grid convention, row/block work ownership
  and existing tree-selection policy.

* Automatic nearest-neighbor density, regression and conditional-density
  recovery now uses each criterion's deleted-sample count domain while
  retaining the full-sample CVAIC endpoint. Explicit and genuinely extended
  starts, selected solvers and healthy search paths are unchanged.

* Generalized-neighbor deleted objectives and native leave-one-out hats now
  decode boundary counts on the deleted donor sample, retaining ordinary-fit
  geometry and existing computational owners. Invalid beta density CVML
  geometry fails closed rather than returning an objective of zero.

* Automatic nearest-neighbor smooth-coefficient bandwidth searches can recover
  a valid deleted-sample boundary when tied conditioning values invalidate the
  initial search points. Explicit starting bandwidths retain their behavior.

* Beta-kernel local-polynomial derivative hats preserve supported external
  rows when another query has all-zero kernel weights, returning undefined
  rows as NA with the existing warning. Complete-operator requests remain
  strict.

* Local-constant conditional-distribution cross-validation uses absolute
  adaptive-neighbor predictor weights consistently in its legacy sums and
  kernel-weight exports. Fixed and generalized-neighbor weight units are
  unchanged. Separate nearest-neighbor fold-geometry limitations remain.

* Private raw delete-one nearest-neighbor kernel geometry no longer consults
  ambient bandwidth-scaling state. Cold beta folded-weight calls use the
  existing context-aware distance owner directly; ordinary public kernel-sum
  conventions and existing computational ownership are unchanged.

* Ordered Li-Racine weights and Racine-Li-Yan kernels accept finite fractional
  numeric distances in their original units. Declared unused levels remain
  part of RLY normalization; alphanumeric levels retain declared ranks.
  Wang-van Ryzin and normalized Li-Racine retain their unit-lattice contract.
  Singular elementwise bandwidth scores at zero smoothing report a precise
  condition. RLY vector evaluation no longer allocates an unused power table.

* Native bandwidth preparation and local-polynomial operators retain declared
  categorical training support, including unused levels. Density,
  distribution, conditional and regression/LSQ routes now use the same
  support metadata as public fits.

* Conditional-quantile inversion accuracy scales with the observed response
  span rather than its absolute location, so translating the response no
  longer loosens the requested accuracy.

* Mixed-data beta copulas retain the applicable kernel components in each
  marginal, including categorical-only marginals, without discarding the
  beta kernels used by continuous variables.

* Density bandwidth documentation now explains how to identify and interpret
  solutions at the fixed-bandwidth search floor, including beta CVLS
  sensitivity to endpoint observations.

* Non-beta first local-polynomial hat derivatives now retain undefined
  external rows as NA with the existing row warning, consistently across
  matrix, application and prediction. Internal finite-operator consumers
  remain strict; native beta partial-row behavior is unchanged.

* Automatic nearest-neighbor distribution-bandwidth recovery now respects
  the deleted-sample count domain. It can recover a valid ordinary candidate
  from a failed boundary start without replacing explicit or extended starts.

* Conditional-mode omission metadata now describes the active evaluation
  result; separate training omission fields retain the training history.
  Estimates, standard errors, class probabilities and NA padding are unchanged.

* Kernel-sum documentation distinguishes zero-diagonal leave-one-out sums
  from deleted-sample nearest-neighbor geometry, and raw weight exports from
  bandwidth-normalized sums. Estimator documentation now explains the existing
  generalized-neighbor distinction between training and explicit evaluation rows.

* Higher and mixed local-polynomial hat derivatives with automatic ridging now
  share the canonical LP solve-admission policy. This avoids spurious ridge
  after translating predictors. Explicit positive ridge remains unchanged;
  default operators report all-zero external rows as unavailable.

* Local polynomial uncertainty recognizes exact interpolation when an accepted
  unregularized local design has exactly as many supported observations as
  coefficients. Unidentified uncertainty is reported as unavailable instead of
  amplifying solve roundoff; point estimates and ridge admission are unchanged.

* Conditional quantile plots preserve fitted numerical extraction controls,
  including resampled fits and bootstrap centers. Explicit plot controls
  take precedence; default and legacy-object behavior is unchanged.

* Smooth-coefficient level hat operators include the fitted model's intercept
  correction after a positive ridge, including categorical-profile and
  local-polynomial paths. Zero-ridge solves and ridge selection are unchanged.

* Partially linear hat operators apply the fitter's identification check to
  the residualized linear design, including subtraction error. Unidentified
  models now fail consistently instead of returning row-order-dependent hats.

* Location-scale quantile fits retain trained formula transformations for
  prediction and scalar children of vector-tau fits. Native exdat is already
  in model coordinates and takes precedence without a second transformation.

* Conditional quantile fits retain their resolved numerical extraction controls.
  Prediction reuses them unless explicitly overridden; old saved objects without
  these controls retain the historical defaults.

* Single-index fitting, prediction and inference carry retained fixed kernel
  endpoints into their regression and kernel-sum consumers, including bootstrap
  resamples. The fitted domain no longer disappears during these calls.

* The Gaussian2 boundary-density helper preserves its selected second-order
  kernel at very large bandwidths instead of switching to a uniform density.
  This removes an artificial discontinuity in fitted values and CV objectives.

* Copula plots align surfaces and uncertainty layers with their probability
  axes when supplied grids are unsorted. Stored fits and returned plot data
  retain their original row order.

* Copula probability-grid evaluation keeps marginal inverse grids inside the
  fitted kernel domain, allowing bounded kernels to retain their valid support
  when the requested grid extension reaches beyond it.

* Mixed-data copula marginals now retain the bounds of their original
  variable, independently of where categorical columns occur in the data.

* Smooth-coefficient evaluation preserves positional native predictor columns
  when preparing large-bandwidth eligibility geometry; harmless vector/frame
  naming differences no longer stop evaluation or bypass qualified shortcuts.

* Regression hat applications use the regression's retained numerical coding
  when borrowing an implicit factor response from a formula, bandwidth object
  or fitted model. Explicit numerical response payloads are unchanged.

* Quantile prediction preserves native exdat coordinates for formula fits,
  avoiding repeated predictor transformations. Native evaluation retains
  precedence over raw formula newdata.

* Regression hat applications now pass resolved degree, basis and Bernstein
  overrides consistently to every numerical owner, including single-response
  shortcuts and retained-operator prediction.

* Prediction from a formula regression hat operator now evaluates newdata
  with the trained predictor transformations. Explicit native exdat takes
  precedence and is not transformed again.

* Recomputing a retained regression hat operator now preserves whether its
  targets are training observations or external queries. Explicit new grids
  remain external even when their coordinates equal the training data.

* npreghat() now dispatches fitted regression objects to their registered
  method before inspecting saved calls, preserving retained training inputs
  and ordinary subclass dispatch.

* Quantile formula refits now consume the resolved missing-data policy before
  filtering numerical controls, allowing the documented replacement-data
  na.action override without weakening native argument validation.

* Explicit gradient orders are now checked consistently by regression fits,
  regression plots and stored gradient/gradient-SE extractors. Unsupported
  higher orders no longer silently return first derivatives; existing LP
  partial availability and categorical first differences are preserved.

* Conditional-mode class support now remains owned by the fitted training
  data. Evaluation response levels contribute to diagnostic tables without
  adding probability columns or changing fitted class probabilities.

* Mixed convolution/tree kernel-weight exports now clear scratch entries
  outside the certified product-kernel support. Sums were already restricted
  to that support; raw matrices now agree with their kernel products.

* Compact generalized nearest-neighbor convolution tree queries now account
  for both kernels' bandwidths, retaining donors whose supports overlap.
  Fixed-bandwidth and adaptive-neighbor computational choices are unchanged.

* Categorical kernel-sum scores now differentiate the requested convolution
  or ordered CDF operator, rather than substituting the ordinary-kernel
  derivative. Li--Racine convolutions also use their retained-support owner
  when scores or categorical profile caching are active.

* Ordered factors with decimal-offset integer-spaced labels now preserve their
  intended distances in kernels, bandwidth scores, CDFs and overlap sums.
  Validation tolerates bounded floating-point representation error while
  continuing to reject genuinely fractional spacing; numeric gaps are retained.

* Conditional-moment test covariance factors now use physical bandwidths,
  including for scaled objects. Quantile-test density weights use normalized
  kernel sums, preserving row/donor-local nearest-neighbor radii rather than
  treating neighbor counts as kernel volumes. This can change non-density-
  weighted NN tests; fixed-bandwidth standardized tests retain their existing
  cancellation of the common scale factor, up to roundoff.

* Copula marginal estimators now use the physical bandwidth retained by the
  joint bandwidth object instead of interpreting scale factors as bandwidths.
  This aligns probability coordinates, density ratios, uncertainty and
  bootstrap marginal operators for scaled and equivalent unscaled objects.

* Unconditional density and distribution fits now consistently use the
  physical bandwidth retained by scaled bandwidth objects. This corrects
  fixed beta-kernel estimates and dependent copula results when
  `bwscaling = TRUE`.

* Type-II significance-test reselection now replaces the physical bandwidth
  and matching scale metadata for the tested predictors together. Previously,
  scaled objects could retain the original smoothing parameters after a search.

* Regression fitting now consumes the physical bandwidth retained by its
  bandwidth object, as numerical hat operators already do. This removes a
  legacy scale reconstruction that could fit a different model when
  `bwscaling = TRUE`; unscaled fitting and bandwidth searches are unchanged.

* Generalized-nearest-neighbor distribution CV now evaluates each bandwidth
  neighborhood on the deleted sample. Empirical query rows retain their
  training identity; explicit grids retain external-query geometry. Canonical
  kernel factors are reused across folds, preserving tree and MPI ownership
  in the corresponding package. Corrected criteria can change selected
  generalized-NN bandwidths.

* Automatic fixed-bandwidth beta density CVLS now requires a positive
  `scale.factor.search.lower` when the complete response sample contains
  repeated observations at a declared endpoint. This conservative search
  restriction leaves the default 0.1 floor, manual evaluation, and
  nearest-neighbor searches unchanged; it does not insert a floor or jitter.

* Beta density CVLS quadrature uses the interior-limit contribution of boundary
  observations instead of assigning positive integration weight to isolated
  endpoint spikes. This applies to unconditional and conditional objectives,
  including distributed quadrature. Ordinary PDF values and CV cross terms
  are unchanged; zero-concentration uniform components remain included.

* Beta generalized-NN unconditional density and distribution fits now preserve
  training-query identity when choosing radii, as other continuous kernels do.
  Explicit external queries keep their external geometry; full-fit standard
  errors use the same corrected contributions. This can change training fits.

* Bandwidth selection now reports activity during baseline-reference
  evaluations without advancing optimizer iteration counts. Native progress
  callbacks preserve interrupts instead of consuming Ctrl-C.

* Adaptive nearest-neighbor distribution CV now evaluates donor radii on each
  deleted sample for beta kernels and external grids as well as empirical
  grids. Canonical CDF factors and occurrence-safe exclusion intervals are
  reused across folds. This can change selected adaptive-NN bandwidths;
  fixed and generalized-NN criteria are unchanged by this repair.

* Ordered Racine--Li--Yan kernels now share one retained-support implementation
  for mass, CDF, bandwidth scores, overlaps and categorical profile plots.
  Regression CV paired moments preserve directed donor weights and row-specific
  AIC diagonals. Near-upper shortcuts certify the full retained support.
  Numeric ordered levels require finite increasing integer spacings; use
  nonnumeric ordinal labels when rank spacing is intended.

* Automatic density-equality test bandwidths now use one pooled-sample search
  and the harmonic-mean sample size to rescale its bandwidth. The same physical
  bandwidth is held fixed in both samples and all bootstrap draws. Categorical
  rescaling is constrained to the kernel's legal maximum, with a warning when
  reached. Supplied bandwidths and manual scale-factor interpretation are
  unchanged; this is a selection policy, not a claim of power optimality.

* Adaptive nearest-neighbor leave-one-out kernel moments no longer risk a crash
  from uninitialized state, including smooth-coefficient bandwidth selection.

* Convolution kernel sums now compute requested derivative and categorical
  replacement products instead of reading uninitialized permutation buffers.
  Raw exported kernel-weight matrices also honor their documented independence
  from sum normalization, including adaptive nearest-neighbor bandwidths.

* Scaled bandwidth objects now supply their retained physical bandwidths to
  categorical profile fits, ridge scaling, hat operators and wild bootstrap
  helpers. Fixed bandwidths equal to one are no longer misclassified as
  nearest-neighbor geometry by the regression validator.

* Conditional properness certification recognizes normalized unordered
  Li--Racine response kernels. Bootstrap surfaces use the same properness
  certificate as the point estimator instead of renormalizing an already
  proper estimator on its finite plotting grid. Degree-zero conditional
  local-polynomial bootstrap helpers share the local-constant operator.

* Quantile and classification formula fits preserve the effective evaluation
  rows when native evaluation arguments override newdata, without applying
  training omissions to the new sample. Their default methods retain positional
  manual bandwidths. Joint native inputs with inconsistent row counts are
  rejected before recycling in quantile, classification, significance-test
  and smooth-coefficient bandwidth routes.

* Scalar bounded conditional CVLS now honors explicit hybrid quadrature ratios;
  uniform components retain the exact declared support endpoints. Gaussian2
  auxiliary boundary-density moments use stable arithmetic at large bandwidths.
  IV derivative initialization rejects contradictory supplied values at repeated
  coordinates before bandwidth preparation. Bootstrap plots distinguish a
  legitimate NA-labelled factor level from a missing observation.

* Uniform-kernel convolution no longer divides by its bandwidths twice.
  This corrects overlap sums and density CVLS objectives at non-unit bandwidths.

* Gaussian eighth-order and Epanechnikov convolution calculations now use
  location-independent polynomial arithmetic. This prevents translation-driven
  cancellation from corrupting kernel overlaps and density CVLS objectives.
  Kernel definitions, bandwidth geometry and computational routing are unchanged.

* Kernel sums with multiple weight and response columns now return the
  documented weight-by-response layout for every kernel and operator.
  Smooth-coefficient backfitting and IV moment consumers use that same layout.

* Kernel-sum bandwidth scores now differentiate the Wang--van Ryzin and
  normalized ordered Li--Racine kernels correctly at all category distances.
  Normal kernel sums and categorical contrast definitions are unchanged.

* Uncompressed density CVLS now squares sample counts in floating-point
  arithmetic, avoiding integer overflow that could corrupt bandwidth selection
  for large samples. Kernel calculations and computational routing are unchanged.

* Bounded conditional-density nearest-neighbour CVLS now evaluates both
  objective terms using the sample with the held-out occurrence removed.
  Response quadrature reuses the existing nodes and bounded work tiles while
  selecting fold-specific radii. This can change selected NN bandwidths;
  fixed bandwidths and unbounded analytic-convolution objectives are unchanged.

* Nearest-neighbour leave-one-out regression beta hats and smooth-coefficient
  fits, hats, uncertainty and CV moments now construct raw weights after
  excluding the held-out occurrence. Adaptive-NN LC moment calculations also
  retain donor-specific bandwidth normalization. Ordinary public kernel-sum
  zero-diagonal semantics and fixed-bandwidth geometry are unchanged.

* Kernel sums with multiple weight columns and an omitted response now retain
  the correct output dimensions, including permutation sums. An omitted
  response acts as a column of ones; kernel calculations are unchanged.

* Beta nearest-neighbour regression CVLS and density CVML use neighbourhoods
  with the held-out occurrence excluded. Unconditional density CVLS applies
  the same correction to its cross term while retaining the full-sample
  integrated squared density. Fixed-bandwidth calculations are unchanged.

* Conditional properization no longer treats LC/degree-zero fitting alone
  as a guarantee of a proper density or CDF. The shortcut also checks kernel
  positivity, response normalization and nearest-neighbour orientation.

* Fixed categorical class-probability effect plots now apply the fitted
  probability projection at both category endpoints before subtraction,
  matching public class-probability predictions.

* Conditional fitted objects retain their resolved properization controls for
  prediction and plotting. Explicit consumer controls take precedence; older
  objects without retained controls keep the historical default behavior.

* Conditional density/distribution properization groups conditioning rows by
  exact identity, avoiding merged slices after large coordinate translations
  or numeric factor-label rounding.

* Categorical recoding and declared-support reconstruction preserve valid
  factor levels labelled NA without converting them to missing observations.
  Actual missing factor codes remain missing.

* Formula fits with native evaluation data no longer restore training
  exclusions into an unrelated evaluation grid. Native evaluation arguments
  consistently take precedence over formula newdata in these routes.

* Smooth-coefficient and partially linear fits honor separately supplied
  evaluation Z data when evaluation X defaults to the training data, and
  validate the paired evaluation row counts.

* Explicit evaluation responses at training locations follow the training-row
  omission mask and native tree ordering. Single-index fits no longer overwrite those supplied
  responses when calculating goodness-of-fit statistics.

* Retained-bandwidth conditional, single-index, partially linear and
  smooth-coefficient fits reject mismatched paired data rows before recycling
  can alter the sample.

* Uniform native quadrature grids retain the exact supplied support endpoints,
  preventing roundoff from dropping a bounded-kernel endpoint contribution in
  conditional and unconditional density cross-validation.

* Mixed categorical/beta density and distribution operators now pass the
  categorical-compression option to their shared native kernel-sum owner.

* Mixed local-polynomial derivatives correctly differentiate every coordinate,
  including constant factors and degree-zero axes. All-zero derivative orders
  return the value basis, including its intercept.

* Conditional fits, operators and bootstrap preparation consistently resolve
  retained scale factors to physical bandwidths, including the separate X/Y
  calculations used by local-polynomial fits. Resampling does not rescale the
  selected smoothing parameters. Bias-corrected bootstrap pilots preserve
  physical/scale-factor metadata in regression, conditional, smooth-coefficient
  and partially linear bandwidth objects.

* Fourth- and sixth-order Gaussian convolution calculations now preserve
  translation and bandwidth-unit invariance. This corrects affected kernel
  sums and density cross-validation calculations.

* Auxiliary beta density kernels consistently use dimensionless bandwidths on
  normalized support. Beta2 rejects h > 1/4 and inadmissible explicit grids;
  automatic starts and optimizer bounds respect that same domain.

* Sampled cumulative quadrature is anchored at zero and invariant to repeated
  evaluation coordinates. Endpoint correction applies to uniform-grid totals
  only; irregular-grid totals use the nonuniform trapezoidal rule. Training
  multiplicities are retained, and the IV derivative iteration reuses prepared
  shared integration geometry. Entropy's corrected-total arithmetic is unchanged.

* Auxiliary boundary and shape-constrained distribution estimates integrate
  the retained fitted density from the true support lower bound, not the first
  requested evaluation point. Infinite tails are not truncated and shape QPs
  are not refitted for integration. Boundary proper=TRUE uses positive-part
  whole-support normalization instead of a constant shift that could create
  infinite mass. Sparse-grid CDF queries can have additional integration cost;
  CDF integration adds no cross-validation work or raw-density refitting.

* The auxiliary Gaussian2 boundary kernel evaluates its removable moment-ratio
  limit at and near the support midpoint, including narrow-kernel ordinate
  underflow, instead of returning NaN. The mathematical kernel is unchanged.

* Conditional-mean and quantile specification bootstrap refits retain model
  observation weights and, for mean models, offsets. Weighted rq models use
  their retained raw model frame rather than weight-scaled x/y components.
  Quantile tests reject a tau inconsistent with the fitted model, and both
  formula tests retain the correct omitted-row indices.

* Auxiliary univariate density/CDF standard errors use the training sample
  size rather than the number of evaluation points. Shape-constrained density
  fits accept one-sided density bounds correctly; integral.equal uses the
  rightmost coordinate regardless of input order and rescales raw derivatives
  consistently with the density (log derivatives are unchanged).

* Density-equality tests freeze one common declared categorical support,
  including unused levels and supplied bandwidth metadata, in every observed
  and bootstrap contraction. Raw numeric bandwidths no longer acquire different
  category normalizations in the two samples. Active ordered Racine-Li-Yan
  kernels are rejected in this test only. Cross-sample normalization no longer
  overflows when the product of integer sample sizes exceeds 2^31-1.

* Conditional-mean and quantile specification tests accept manual bandwidths
  without duplicate kernel arguments, retain selected kernel metadata in all
  contractions, and keep na.exclude models on their compact estimation sample.
  They validate actual model components rather than the spelling of the model
  call; incompatible test/model response samples are rejected explicitly.

* Formula column extraction and retained-role replay support quoted predictor
  names without silently renaming the data. Pipe-formula splitting distinguishes
  actual operators from punctuation inside quoted names or transformations.
  Fitting, refitting, prediction, significance tests and plotting share the same
  symbolic-name interpretation; formulas are not re-evaluated to decode names.

* Density-equality tests now use one common bandwidth in both samples and
  all bootstrap replications. Automatic selection uses pooled LSCV with the
  harmonic-mean reference-size adjustment described above; explicit CV-method
  choices remain available. One supplied bandwidth is reused without a search;
  incompatible two-bandwidth inputs are rejected.
  This test now explicitly rejects nearest-neighbour bandwidths and
  boundary-normalized/beta kernels, whose directed cross sums need not be
  symmetric as its statistic and studentizer require. Fixed unnormalized
  Gaussian, Epanechnikov and uniform kernels remain supported. No other
  estimator's kernel or nearest-neighbour capabilities change.

* Bootstrap upper-tail p-values count ties as at least as extreme in all
  inference routines, consistently with significance tests. They remain
  ordinary bootstrap proportions (no add-one correction); zero p-values
  remain possible. Discrete/categorical tied cases may now report larger
  p-values.

* Smooth-coefficient hat application with categorical profile compression
  now expresses a positive ridge in the same kernel-weight units as the
  matrix output. Bandwidth selection and zero-ridge arithmetic are unchanged.

* Native data-first fitting calls such as `npudens(x)` again pass observations
  to bandwidth selection instead of interpreting them as bandwidths. Explicit
  manual-bandwidth and named training-data calls retain their interpretation.
  Long univariate input expressions no longer fail while constructing a
  single variable label.

* Univariate entropy equality and symmetry tests retain declared categorical
  support, including unused levels, throughout probability summation and
  bootstrap resampling. Equality tests use common category metadata; ambiguous
  qualitative ordered supports require an explicit common level declaration.
  Corrected categorical statistics, selected bandwidths and p-values may differ.
  A shared plug-in selector now evaluates each sample's own probabilities,
  fixing equality tests that supply only the first bandwidth. Independent
  missing-value omission is completed before overlap checks, and obsolete
  omission attributes no longer break numeric entropy integration.
  Categorical symmetry now rejects an observed or bootstrap reflection outside
  its declared support instead of silently losing those observations.

* Bootstrap progress retains the latest completed work between heartbeats.
  Operator rows and blocks are distinguished from completed replications;
  conditional-moment tests credit evaluated statistics, not generated draws.

* Serial-dependence tests reject missing time points rather than silently
  joining observations across gaps. Supply a contiguous complete segment
  explicitly when the original series contains missing values.

* Kernel sums, hat operators and quantile evaluation now distinguish current
  missing rows from historical na.action attributes on already-clean data.
  Supplying na.omit/na.exclude output no longer deletes unrelated rows a
  second time. Formula-result omission restoration remains unchanged.

* Density-equality tests retain an individually supplied bandwidth while
  selecting only the missing counterpart. Incomplete rows are omitted once
  per independent sample before bandwidth selection and bootstrap resampling,
  avoiding inconsistent kernel-sum/count dimensions; at least two complete
  observations are required in each sample.

* Significance tests keep their progress heartbeat active during long native
  bootstrap tiles and reuse an already-computed unrestricted fit when preparing
  streamed residuals. Bootstrap draws and test definitions are unchanged.

* Density-equality tests display activity while computing the observed statistic
  after the bootstrap replications finish.

* Reusing a single-index bandwidth object created by a one-call formula fit
  now preserves its training omission map. With na.exclude, fitted values,
  standard errors, gradients and residuals retain their original row positions,
  including residuals requested after fitting.

* Smooth-coefficient prediction and location-scale quantile SE extraction now
  preserve a computed NA standard error for a single unsupported evaluation
  row, as they already did for larger grids. Unrequested uncertainty still
  gives a no-search refit message.

* Location-scale quantile prediction and retained-bandwidth refits now restore
  omitted rows consistently, including standard errors and multi-quantile
  residuals. Explicit replacement samples keep their own row positions.
  Scalar numeric I(...) terms also work in location-scale quantile evaluation
  and kernel-sum formula evaluation, as they do when fitting.

* Formula regression, single-index and smooth-coefficient residuals now retain
  training-row ownership when an external evaluation grid has different missing
  rows. Smooth-coefficient coefficients and effects are restored alongside
  fitted values under na.exclude. External smooth-coefficient residual requests
  no longer depend on requesting standard errors.

* Smooth-coefficient prediction and lazy residual extraction retain the fit's
  iteration, tolerance, iteration limit and leave-one-out controls. Explicit
  prediction controls still take precedence.

* Scalar numeric I(...) formula terms now use the ordinary numeric predictor
  type, without admitting matrix-valued terms or changing the stored formula.

* Smooth-coefficient wild-bootstrap coefficient intervals now reuse the same
  categorical kernel normalization as their fitted moment equations. Density
  profile bootstrap intervals also retain the estimator's normalized
  Li--Racine kernels for both unordered and ordered factors. Bootstrap draws,
  ordinary fitted coefficients and bandwidth selection are unchanged.

* Local-polynomial fitting, operator and bootstrap rank certificates now count
  exact repeated design rows once, rather than treating duplicate observations
  as independent directions. Sparse supported designs use the existing ridge
  policy consistently across fits, hats, inference and counted bootstrap draws.
  Identity preparation is shared within each fitting invocation; bandwidth
  search and numerical ridge thresholds are unchanged.

* Significance-test progress now shows the active predictor's position and
  completed bootstrap replications, with call-wide elapsed time and ETA.
  Analytic non-rejection is labelled as skipped bootstrap work.

* Partial response replacement with formula-derived bandwidth objects now uses
  retained-sample row order, before jointly applying the effective NA policy.
  Refits retain the replacement sample without changing the original bandwidth
  object. Original-data replacement belongs in `data`; it cannot be combined
  with a partial response. Bandwidth selection is not repeated.

* Native-object standard `newdata` now matches named columns to the retained
  variable names in fitting and prediction. Reordered or extra columns are
  handled consistently; missing, duplicate or ambiguous names are rejected.
  Explicit native evaluation arguments remain positional and take precedence.
  Unnamed single-role matrix/vector inputs remain positional. Copula
  `newdata` uses the documented probability-coordinate names.

* `npreghat` rejects explicit nonzero `ridge` requests on computational
  owners that cannot honor them. Default automatic numerical ridging and
  supported generic mixed-derivative ridge requests are unchanged.

* Density/distribution hat and bootstrap adapters now preserve normalized
  Li--Racine ordered kernels, including conditional response kernels.
  This corrects affected ordered-factor operators, bootstrap intervals and
  conditional local-polynomial contrast/higher-derivative uncertainty.
  Raw regression and kernel-sum ordered-kernel conventions are unchanged.

* Hat derivative selectors now validate orders before integer conversion
  and honor named coordinates and explicit prediction overrides. Index hats
  accept matrix responses consistently and retain singleton matrix dimensions.
  Supported explicit hat ridges are retained when rebuilding an operator.

* Fits reusing retained native bandwidth data now give explicit replacement responses
  precedence over retained defaults. Copula predictions preserve inversion
  settings; conditional-mode predictions preserve the requested probability
  projection. Copula joint and marginal calculations use a common named column
  order and the actual replacement sample size.

* Formula kernel-sum weights and location-scale regression scales now follow
  the response/predictor subset, time alignment and missing-value selection.
  Full-length auxiliary inputs refer to the original sample; unambiguously
  preselected shorter inputs remain supported.

* Native one-call fits pass their already evaluated training data to the
  bandwidth constructor. Side-effecting or random data expressions are no
  longer evaluated again for that handoff; fitting and retained bandwidth
  data use the same realization, including MPI dispatch.

* Local-polynomial fitting, hats and related inference now consistently
  apply the existing ridge rule when fewer nonzero-weight donor rows
  contribute than there are basis coefficients. This corrects unstable
  fit/hat disagreement at structurally underdetermined rows, including
  mean-only fits. Bootstrap rank counts use independent donor rows rather
  than resampling multiplicities. Bandwidth-search objectives and the ridge
  magnitude and numerical thresholds are unchanged.

* Conditional gradient bootstrap preparation reuses the canonical internal
  bandwidth builder and finalizer instead of repeating public constructor setup.
  Each resample still uses its own geometry and the same numerical fit.

* Bootstrap plot summaries now align categorical draw columns with their
  actual evaluation labels, including declared but unobserved factor levels.
  This fixes unused-level plots across regression, conditional and
  unconditional, partially linear and smooth-coefficient families without
  dropping levels or changing kernels, draws or intervals.

* Native bandwidth-object fits now honor `newdata` using the corresponding
  prediction interface's evaluation roles instead of silently returning the
  training fit. Conditional and semiparametric frames must identify their
  roles by retained column names; explicit native evaluation arguments keep
  precedence.

* Exact wild regression bootstrap helpers restrict LC derivative-operator
  reuse to its intended continuous LC/LP0 first-derivative requests, and align
  missing evaluation rows before blocking. Other derivative requests retain
  their per-response estimator, including the existing beta-kernel path.

* Ordered kernel-sum contrasts validate evaluation category indices before
  native work and initialize the cached lookup on the first row. Invalid
  contrast categories now give a clear error instead of an invalid memory
  read; ordinary ordered-kernel extrapolation is unchanged.

* Adaptive nearest-neighbor leave-one-out regression hats now supply matching
  single-evaluation-point tree geometry to compact-kernel traversal, avoiding
  an out-of-bounds read without disabling tree evaluation or changing radii.

* Interval cap widths now use horizontal display spacing consistently for
  level, gradient and partially-linear coefficient plots, independent of the
  response-axis range. Exact-zero intervals, interval endpoints and bootstrap
  draws are unchanged.

* Smooth-coefficient wild-bootstrap coefficient-gradient plots now reuse a
  direct, blocked coefficient projection across bootstrap responses. The
  projection retains the fitter's normalized moments, accepted ridge and
  intercept correction for LC/LL/LP and fixed/generalized/adaptive NN
  bandwidths. Pilots, draws, RNG order and intervals are preserved to numerical
  rounding. Mean plots and pairs/block resampling retain their existing owners.

* Local-constant and degree-zero local-polynomial continuous-gradient wild
  bootstrap plots with Gaussian, Epanechnikov and uniform kernels now reuse
  the exact derivative operator across responses,
  including local-smoothing quantile plots. Operator rows are blocked for
  large grids. Bootstrap draws, pilots, intervals and RNG ordering are unchanged
  to numerical rounding; pairs/block resampling retains its own geometry.
  Beta kernels retain their per-response calculation because endpoint
  derivative cancellation can depend on the response.
  The canonical LC derivative operator now returns NA at external rows whose
  computed kernel weights are all zero, matching direct fitting; complete-row
  internal consumers still reject them. Invalid radii and other nonfinite
  results are not converted into partial results.

* Plot interval ranges now align NA-padded categorical panels with their actual
  interval rows. Unequal partially linear panel sizes no longer generate a
  recycling warning with `band="all"`; interval values and axis ranges are
  unchanged. Inconsistent finite panel data are rejected rather than recycled.

* Formula training now resolves character NA policies through the same owner
  as `stats::model.frame`, so an unrelated caller-local `na.omit` no longer
  changes the sample. Explicit NA function objects retain their own behavior.
  With replacement training data, an explicitly supplied `na.action` takes
  precedence over the saved policy, including `NULL`, and the fitted object
  retains the effective policy. Native complete-case exclusions remaining
  after a permissive policy are now recorded with the retained sample rather
  than lost from formula omission bookkeeping.

* Internal native bandwidth children without an original construction call no
  longer acquire an executable placeholder call. Their retained-sample refits,
  predictions and plots remain supported; `update()` on a call-less child now
  fails instead of potentially using unrelated caller data and selecting new
  bandwidths.

* Fits with replacement training data now retain that fitted sample for later
  refits, prediction, plotting and inference. Native bandwidth constructors for
  regression, density/distribution and semiparametric methods also retain their
  training values, so rebinding caller variables does not change later results.
  Formula and native sample ownership no longer compete in LSQ scale-pilot
  children. Objects saved without retained samples still need their original
  data or explicit replacement inputs.

* Positional numeric bandwidths in one-call fitting are now honored like named
  bandwidths, without unintended bandwidth selection. Existing inline formula
  and wrapper syntax is unchanged.

* Ordered categorical contrasts correctly initialize the category lookup when
  the first evaluation category is zero. This corrects affected local-constant
  gradients, standard errors and their downstream bootstrap/test consumers.

* Interval-cap widths now use display geometry consistently across response
  scales, including singleton evaluation points.
  Partially linear independent-scale plots with data overlays now use the
  correct panel-local fitted values and bootstrap summaries.

* Partially linear and local-smoothing quantile child bandwidths retain their
  own native samples and responses for later use. This is sample ownership,
  not a general environment-compaction guarantee: user formula and custom
  function environments remain available and can retain unrelated caller
  objects in saved files. Retained samples also increase object storage;
  derived children may retain separate copies. Serialized size can exceed
  the size reported by `object.size()`.

* Formula bandwidth objects now retain their prepared training sample and NA
  policy. Rebinding constructor variables (for example in a model-building
  loop) no longer silently changes later refits, predictions, plots or tests.
  Explicit training-data overrides remain available. Objects created by older
  versions without retained samples still require their original data bindings.
  Single-index bandwidth plotting uses the same saved-data owner.

* Arguments forwarded through `...` wrappers are honored as in direct calls:
  an explicit `bws` is no longer ignored in favor of bandwidth search.
  Saved formula objects also retain forwarded or caller-local `na.action`
  settings when refitting, predicting, plotting or running significance tests.

* Formula constructors supplied inline are evaluated once per invocation,
  including the one-call bandwidth/fit route. Existing formula values and
  literal formulas remain supported; callers need not construct a formula
  in a separate statement. Forwarded calls retain lazy subset evaluation.

* Streamed non-studentized categorical significance-test bootstrap statistics
  now use the same generalized-nearest-neighbour training-row geometry as
  the observed statistic and public regression refits. This can change
  affected bootstrap statistics and p-values; pivot=TRUE is unchanged.

* Conditional density/distribution formulas reject unsupported transformed
  responses instead of silently discarding their transformations. Prepare the
  transformed response in the data, or use the native data interface. Formula
  routes also reject unsupported RHS offset specials before evaluation;
  supported response/predictor transformations and ordinary variable names
  remain unchanged.

* Automatic regression wild-bootstrap pilots now use fitted training-row
  geometry whether plotting starts from a fitted model or its bandwidth
  object. This corrects generalized-nearest-neighbour bands that previously
  depended on that entry route. Explicitly supplied pilots and external
  evaluation geometry are unchanged.

* The significance-test documentation now notes that studentization can
  improve size behavior without uniformly improving finite-sample calibration
  or power. The default remains pivot=TRUE.

* Formula fitting defers data-dependent subset expressions to the model-frame
  owner rather than forcing them during argument dispatch. One-call fits and
  constructor-first fits now use the same selected rows.

* Partially linear prediction standard errors use the prediction's actual
  retained training rows, including data overrides, and preserve na.exclude
  padding of evaluation rows. Public residual padding is unchanged.

* Vector-quantile least-squares predictions prepare formula newdata once and
  share that same evaluation sample across quantiles.

* Partially linear formulas now prepare their response, linear regressors and
  smoothing variables together once, retaining trained transforms and a common
  subset/time-aligned sample through fitting, prediction and plot preparation.

* Regression hat/significance helpers and smooth-coefficient plots now resolve
  saved training-data expressions in their original call owner, including
  bandwidth objects constructed inside wrappers.

* Copula formula construction and fitting now share the same prepared sample,
  including transformed, subsetted and time-aligned inputs. Automatic grid
  dimension validation still precedes bandwidth search.

* Smooth-coefficient formula construction, fitting and evaluation reuse a
  single prepared sample per transaction and retain portable trained terms.

* Single-index formula construction and prediction now use the shared
  single-evaluation preparation route and portable trained prediction terms.
  Named formula calls no longer depend on argument order or fail after
  bandwidth construction.

* Quantile and conditional-mode one-call formula fits also reuse their prepared
  training sample. Conditional-mode response validation still occurs before
  bandwidth selection. Saved-object refits now reuse the retained training
  sample; explicit evaluation data still use the trained formula terms.

* Conditional density/distribution formula fitting now reuses its prepared
  training sample, avoiding repeated evaluation of predictor expressions.
  Saved conditional formula terms retain portable prediction metadata for
  refitting, quantile/mode evaluation and plotting.

* Fixed/geometric plot-bootstrap draws now use a replicate-ordered random
  stream, independent of processing chunk size and MPI worker count. The
  block law is unchanged, but historical seeded bands may change. The stream
  follows the previous single-replicate ordering; IID/wild draws are unchanged.

* Generalized-NN fitted-row categorical effects and their standard errors now
  retain training identity for local-linear and local-polynomial regression,
  as local-constant regression already does. Explicit evaluation effects keep
  external-query radii. Fitted means and continuous derivatives are unchanged.

* Time-series formula alignment preserves the columns belonging to each
  matrix-valued expression, preventing a later formula variable from being
  replaced by a column of an earlier expression. Existing restrictions on
  matrix predictors are unchanged.

* Regression and unconditional density/distribution formula fits evaluate
  each training expression once per constructor/fit operation. This avoids
  inconsistent samples from random or stateful expressions and retains
  prediction-transform metadata. Later refits still evaluate their inputs.
  Stored formula objects for these families can be used by either package
  without loading the package that created them.
  Unconditional refits, predictions and plots resolve function-local data
  through the retained bandwidth call, including after MPI fitting.
  Automatic regression calls using `formula =` follow the same preparation
  and fitting path as positional formulas.

* Qualified Apple Silicon builds accelerate compensated local-polynomial
  regression influence calculations used by requested standard errors and
  studentized tests. The calculation, defaults and bootstrap random-number
  sequence are unchanged; other builds retain the existing scalar path.

* Copula standard errors and plots recover physical evaluation coordinates
  independently of sanitized output column names. Transformed, lagged and
  reserved predictor names no longer cause plotting failures or substitute
  training/probability coordinates for a requested grid. Public output columns
  and the underlying estimation and uncertainty formulas are unchanged.

* One-call `npreg()` formulas again align lagged time-series variables before
  bandwidth selection, matching `npregbw()` followed by `npreg(bws = ...)`.
  This repairs a 0.70-1 regression that could silently fit unaligned lag
  columns, or reject them as collinear, and affect subsequent inference and
  plots. LC, LL and LP use the same restored formula-constructor route.

* Lagged time-series formula alignment is also applied to single-index,
  least-squares quantile, kernel-sum, copula and IV formula preparation.
  Quantile, conditional-mode and regression-hat evaluation reconstruct the
  aligned predictors from the supplied new data, including transformed
  terms, without requiring columns named after the expressions. Ordinary
  data and explicit native-array alignment conventions are unchanged.

* Ordinary regression inference now normalizes training fitted residuals by
  the norms of their actual residual-smoother rows before forming level,
  derivative and paired categorical-contrast sandwich errors. This can also
  change continuous and categorical studentized test statistics and
  P-values, not just reported standard errors. It uses
  the pre-smoothing residual correction of Ruppert et al. (1997), not a
  separate variance-function smoothing estimator. Requested errors that
  genuinely lack donor residual information are reported as unavailable,
  separately from certified zero influences. No variance floor, row deletion
  or bootstrap-draw filtering is introduced. The correction does not remove
  smoothing bias or guarantee finite-sample calibration under arbitrary
  heteroskedasticity or model selection.

* Requested categorical regression effects use stable paired-influence
  arithmetic with and without standard errors. This can change previously
  cancellation-dominated gradients and downstream test statistics.
  Mean-only fitting and bandwidth search are unchanged. Residual preparation
  is confined to requested inference, but explicit categorical-gradient
  requests also incur more careful arithmetic when SEs are off. The
  correctness changes can increase these requested computations.
  Single-index mean/gradient errors and transformed-response least-squares
  quantile errors inherit the shared correction; separate coefficient
  covariance formulas do not. Help distinguishes this residual correction
  from numerical ridging motivated by Seifert and Gasser (2000).

* Studentized joint `npsigtest()` calls reuse each bootstrap regression fit
  across tested categorical predictors within the existing bounded response
  tile. Statistics, resamples and random-number consumption are unchanged.
  Individual studentized categorical bootstrap fits use the same native
  contrast and response-specific sandwich calculations but omit unrequested
  gradient and contrast-SE consumers. Default full regression output and
  unstandardized tests retain their existing computational paths.

* `npsigtest()` now defaults to `pivot = TRUE` for continuous and
  categorical predictors, including joint tests. Explicit FALSE retains
  unstandardized statistics; NULL is no longer an input mode. Categorical
  studentization uses response-specific paired-contrast sandwich errors,
  with additional computation confined to requested inference. An entirely
  zero observed effect gives analytic non-rejection without bootstrap fits;
  skipped columns and executed counts are explicit. Unexplained undefined
  standard errors produce a diagnostic, never a variance floor or fallback.

* The getting-started guide demonstrates prediction and no-search refitting
  for uncertainty and gradients, distinguishes computation from extraction,
  and summarizes migration from 0.70-5. Package help clarifies practical
  computational limits. No estimator behavior or defaults change.

* Entropy tests refresh elapsed-time progress during native summation and
  quadrature using the shared transient owner. Original-statistic, lag and
  bootstrap phases clear on completion or failure; bootstrap counts advance
  after their statistics finish, without changing resamples or calculations.

* Progress displays use consistent work counts, elapsed-time and ETA wording.
  Narrow consoles preserve complete counters before optional detail, and
  wide-character labels are shortened by their displayed width.

* Requested adaptive-NN standard errors and statistical-test bandwidth phases
  use the existing transient progress owners. Supplied bandwidths no longer
  produce a misleading standalone bandwidth-computation notice.

* Explicit fixed-bandwidth conditional bootstrap bias-corrected plot centers
  now use one mixed-data donor-mixture pilot and its matching PDF/CDF or
  derivative/contrast reference. Categorical conditioning pilot lambdas
  generate normalized transitions on the declared levels; paired contrasts
  share each resample. Level refits now retain the selected LC/LL/LP estimator
  instead of substituting an LC ratio for LL/LP. These corrections can change
  explicitly requested bias centers, including continuous-only panels.
  Ordinary fitted values, estimate-centered variability draws and search are
  unchanged. NN, proper-projection and quantile bias-center exclusions remain.

* Smooth-bootstrap Epanechnikov perturbations now use the same standardized
  kernel scale as estimation. This corrects explicitly requested unconditional
  and conditional density/distribution bias centers; Gaussian/uniform draws
  and RNG consumption are unchanged. It does not change ordinary bootstrap
  variability draws.

* Uniform nearest-neighbour kernels now apply their documented strict support
  boundary to the unstandardized distance. This corrects adaptive- and
  generalized-NN results that could depend on evaluation batching, tree
  partitioning or affine rescaling. Dense and tree consumers share the same
  membership decision, including normal factors in gradient and integral
  products. Fixed-bandwidth arithmetic, other kernels and tree eligibility
  are unchanged; affected uniform-NN fitted values and search objectives may
  change where reciprocal rounding previously admitted a boundary point.

* Unconditional adaptive-nearest-neighbour standard errors use correctly
  aligned transient workspaces on platforms where extended precision requires
  stricter alignment. The formula, precision and SE-off computation are unchanged.

* Local-linear and positive-degree local-polynomial conditional density,
  distribution and quantile fits now provide requested categorical-contrast
  standard errors. Paired endpoint influences retain their covariance;
  quantiles use each endpoint's own quantile and density. This includes
  fixed, generalized-NN and adaptive-NN bandwidths conditional on the realized
  geometry. Existing point estimates, continuous errors and SE-off computation
  are unchanged. Nonsmooth or nonfinite influences remain explicitly unavailable.

* Bootstrap plots using the shared estimator dispatcher
  now retain one display owner and elapsed clock across their scheduled targets,
  with an approximate overall ETA. Copula bootstrap interval preparation uses
  the same scope. Numerical work and standalone-helper behavior are unchanged;
  conditional-mode and LSQ-quantile displays and the frozen bootstrap shortcut retain their
  existing behavior.

* Conditional density/distribution bootstrap gradient plots now retain the
  requested physical predictor column. Categorical panels previously could
  receive a continuous derivative's bootstrap draws, or fail when predictors
  were interleaved. Their intervals now target the displayed first difference.
  Explicit fixed-bandwidth bias centers now use the mixed-data pilot described
  above, including categorical targets.

* Conditional bootstrap gradients reuse the existing evaluator with a private
  single-target demand, avoiding unrequested R-side categorical differences
  and higher-derivative restoration. Native computation, exact NN resampling,
  random draws, public fitting and asymptotic-error behavior are unchanged.
  The mixed-data smooth-bootstrap pilot also reuses this single-target
  evaluator; its separate probability-law corrections are described above.

* Conditional-quantile bootstrap gradients likewise compute only the requested
  R-side categorical quantile contrast. The original selected-CDF inversion,
  density evaluation, native derivatives and requested endpoint inversions are
  retained; full-output fits, tau layouts and LSQ quantile owners are unchanged.

* Conditional density/distribution fitting now consistently retains explicit
  external rows with all computed explanatory-kernel weights zero as NA, with
  one notice. The policy covers LC/LL/LP, requested inference and categorical
  counterfactuals without changing supported estimates or valid zero response
  contributions. Proper projection retains unsupported fixed-X slices while
  projecting supported slices; prediction and plot consumers preserve status.

* Conditional quantile inversion now preserves explicit base-support status
  for external rows across scalar/vector tau, prediction, inference and MPI
  tau blocks. Unsupported outputs are NA with one notice; genuine numerical
  and required internal-evaluation failures remain errors. The point-only
  constant-response shortcut is unchanged and does not preflight X support.
  Conditional-mode evaluation likewise preserves one final notice without
  changing supported class/probability calculations.

* Requested unconditional adaptive-NN density/distribution standard errors now
  include the same-sample influence of donor-specific radii, with a two-sided
  NN rank-spacing pilot. The joint influence retains cross-coordinate and
  categorical-kernel covariance. Gaussian/Epanechnikov density and distribution
  kernels, and uniform kernels, are supported on regular interior
  NN ranks without finite kernel bounds. Uniform density errors include the
  joint marked moving-support faces using sample-quantile/spacing pilots;
  strictly absent faces have zero motion. Unqualified boundary-pilot, bounded,
  extended/saturated-rank or degenerate-spacing cases retain point estimates
  and return NA standard errors with one availability notice. This changes
  previously incomplete ANN uncertainty, not point estimates, bandwidth search
  or SE-off work. Bandwidth selection uncertainty is not included.

* Requested fixed-bandwidth copula density standard errors now include the
  same-sample covariance of the joint and all marginal density estimates.
  The centered ratio influence is reduced in bounded query tiles; MPI keeps
  marginal kernel weights on their owning rank and communicates only SEs.
  Point estimates, bandwidth search, SE-off work, distribution copulas and
  NN uncertainty are unchanged. This correction holds marginal evaluation
  coordinates fixed; it does not include quasi-inverse or bandwidth uncertainty.

* Requested conditional density/distribution LP standard errors now follow
  the fitted regularized operator when explanatory-kernel weights are
  constant. This corrects level and native first-derivative uncertainty
  after ridging; fitted values and point gradients are unchanged.

* Requested scalar conditional adaptive-NN standard errors now use the
  complete empirical ratio influence for levels and existing continuous
  first derivatives, including continuous/mixed responses and categorical-only
  predictors. This corrects their physical bandwidth scaling while holding
  the selected model, predictor geometry and realized response radii fixed;
  sampling or selection uncertainty in these quantities is not included.
  Point estimates and gradients are unchanged. Local-polynomial, beta-kernel
  and paired categorical-contrast errors retain their separate targets.

* Fixed/generalized-NN conditional density and distribution standard errors
  with categorical variables now retain their joint categorical second
  moments. Conditional CDFs and densities without continuous responses also
  retain numerator--denominator covariance, so a constant response-kernel
  contribution has zero variance. Pure-category profile and ordinary paths
  agree. Point values, point gradients, likelihoods, SE-off work and existing
  beta/positive-degree LP/paired categorical-contrast owners are unchanged.
  Broader adaptive-NN variance limitations remain deferred.

* Fixed/generalized-NN unconditional density and distribution standard errors
  with mixed data now retain squared joint categorical factors within their
  existing leading continuous-bandwidth approximation. Constant categorical
  factors consequently scale uncertainty consistently with the estimate.
  Point estimates, likelihoods, SE-off work, pure-continuous/pure-categorical
  conventions, beta kernels and adaptive-NN errors are unchanged.

* Partial-linear fitting and plot/bootstrap helpers now consistently reject
  numerically unidentified residualized linear regressors across BLAS
  implementations. The shared check accounts for column scale and residual
  formation precision while retaining the existing QR, accepted coefficient
  solve and covariance calculation.

* Public single-index fitting and prediction now return `NA` with one notice
  for external rows whose computed kernel weights are all zero, consistently
  across point-only, asymptotic and bootstrap requests. A bootstrap draw that
  loses support makes the affected uncertainty unavailable without replacing
  a supported original point estimate. Required training ratios remain strict;
  supported values, resampling and computational owners are unchanged.

* Requested higher-order conditional density/distribution LP derivative
  standard errors now use the existing local working variance and the actual
  derivative weighting row, including its accepted ridge. Fixed and adaptive-NN
  point owners are preserved; generalized-NN higher-point restrictions remain.
  These errors condition on the selected bandwidths, radii and local map.
  No supplementary uncertainty work is requested with `se = FALSE`.

* Missing stored gradients or standard errors now consistently explain how to
  refit with the retained bandwidth object without repeating bandwidth search.
  Unsupported uncertainty targets retain their separate diagnostic messages.

* Model-specification testing now uses the shared bandwidth progress owner,
  clearing its transient line at completion instead of leaving a permanent
  "Computing bandwidths" message. Test calculations are unchanged.

* Requested conditional local-polynomial standard errors in the general
  fitting path now use the accepted fit's intercept ridge correction.
  Failed covariance blocks or invalid variance calculations return `NA`
  only for affected uncertainty outputs, including beta-response level errors.
  Point estimates, search, `se=FALSE` work and derivative availability
  are unchanged.

* Requested smooth-coefficient standard errors now use the same accepted
  moment system and intercept ridge correction as the returned point fit.
  Supported signed kernels no longer receive a separate covariance-only
  ridge. Fitted values, coefficient effects, search and `se=FALSE` work are
  unchanged; unavailable covariance retains its existing `NA` output policy.

* Copula fitting accepts named `se=FALSE/TRUE`, defaulting to FALSE.
  Requested asymptotic errors are stored at fitting time; `se()` only reads
  stored results and explains how to refit with the bandwidth object when
  errors are absent. Explicit prediction/plot inference requests retain their
  existing uncertainty target; point estimates and NN geometry are unchanged.

* Unconditional density and distribution fits now default to `se=FALSE`,
  skipping uncertainty-only computation while retaining the same point
  estimates. Use `se=TRUE` to retain asymptotic standard errors, or
  `predict(..., se.fit=TRUE)` for prediction errors. `se()` only extracts
  stored results and explains how to reuse the bandwidth object for a refit
  without repeating bandwidth selection. Explicit asymptotic plot and copula
  inference continue to request their required errors.

* Partially linear fitting now defaults to `se=FALSE`, avoiding coefficient
  covariance work unless requested. Use `npplreg(..., se=TRUE)` for
  `coef(..., se=TRUE)` and `vcov()`; missing-output messages show how to
  reuse the bandwidth object without repeating search. Explicit prediction
  standard errors and asymptotic plots continue to request their needed work.

* Conditional density and distribution fitting now have a named-only
  `se=FALSE` default. Disabled errors skip uncertainty-only native work and
  are stored as `NULL`; fitted values and requested gradients are unchanged.
  Use `se=TRUE` for available asymptotic errors, or refit a retained bandwidth
  object to obtain them without repeating search. Prediction with
  `se.fit=TRUE` requests the needed work explicitly.

* Quantile and conditional-mode fits now accept named-only `se=FALSE`.
  Uncertainty-only evaluation and storage are skipped by default, independently
  of the existing gradient controls. Explicit prediction or asymptotic plot
  requests enable their required errors; stored-result extractors give a
  bandwidth-reusing refit instruction when uncertainty was omitted. Point
  quantiles, class selection, probabilities and requested gradients are unchanged.

* Ordinary scalar conditional categorical-gradient standard errors now use
  complete paired same-observation influences, including endpoint covariance.
  Requested errors use an additional streamed inference pass and can increase
  fit time when categorical uncertainty is retained with `gradients=TRUE`.
  Point estimates, gradients, bandwidth selection and computational selectors
  are unchanged. The default `gradients=FALSE` and internal point-only and
  bootstrap consumers do not request this additional work.

* Conditional local-polynomial fits now retain available first-derivative
  standard errors when other requested derivatives exceed their fitted
  degrees. Existing point gradients, higher-order availability and computational
  owners are unchanged; point-only internal consumers request no extra errors.

* Conditional generalized-nearest-neighbor training fits with ordinary
  explanatory kernels and a beta response kernel now exclude the identified
  training occurrence when selecting explanatory radii, retaining it in the
  estimation sum. Affected fitted values, derivatives and standard errors use
  the corrected geometry; explicit evaluation queries and beta explanatory
  kernels retain their existing behavior.

* Partially linear coefficient standard errors and covariance now allow
  conditional heteroskedasticity through a residual-score sandwich, retaining
  the existing finite-sample multiplier. Coefficients, fitted values, bandwidth
  selection and the separate fitted-value prediction-error calculation are
  unchanged. Coefficient intervals and dependent asymptotic plot terms use the
  corrected covariance.

* Scalar conditional categorical-gradient standard errors involving beta
  kernels now use paired same-observation influences at the two endpoints,
  including their covariance. Self-contrasts have zero error by construction;
  point effects, level errors and positive-degree LP covariance are unchanged.

* Scalar conditional beta-kernel influence standard errors now use the global
  sample-covariance factor n/(n-1), correcting an extra division by n in level
  and supported derivative variances. Conditional point estimates, ordinary
  regression and positive-degree local-polynomial covariance are unchanged.

* Scalar conditional beta-response derivative standard errors now allow a
  legitimate zero derivative of the explanatory weight sum. Zero endpoint
  denominators remain invalid; estimates and existing valid derivatives are
  unchanged.

* Scalar conditional Gaussian and Epanechnikov gradient standard errors with
  fixed/generalized-nearest-neighbour bandwidths now use the squared analytic
  derivative moment, not a finite-shift convolution constant. This changes
  uncertainty only; uniform derivative inference and broader adaptive-NN
  variance questions are not resolved by this correction.

* Scalar conditional density and distribution standard errors now use the
  squared-kernel moments of the chosen explanatory kernels and, for density,
  response kernels separately. Mixed kernel families/orders no longer inherit
  the response-kernel moment in every dimension. Point estimates are unchanged.

* Pure-categorical density standard errors now use the empirical variance
  of the complete joint-category kernel contributions, instead of the
  continuous-density approximation. Constant contributions have zero
  standard error, and zero-smoothing probabilities retain the existing
  binomial convention. Density estimates, search, categorical shortcuts,
  and mixed/continuous density uncertainty are unchanged. The correction
  reuses existing kernel rows with bounded moment storage.

* Pure-ordered distribution standard errors now also use the empirical
  variance of complete joint CDF contributions. Finite smoothing no longer
  treats each contribution as a Bernoulli indicator; the zero-smoothing
  convention, distribution estimates, and continuous-variable formulas
  are unchanged.

* Wang–van Ryzin ordered distribution estimates now use the correct cumulative
  probability above each training category. This corrects potentially decreasing
  distribution estimates and affected distribution cross-validation criteria.

* Single-index mean hats and affected bootstrap ratio helpers now normalize
  finite signed and positive-small weight sums directly instead of flooring
  them at machine epsilon. Empty external LC/degree-zero index and adaptive
  conditional-hat rows return NA with one notice; undefined required ratios
  fail without discarding or replacing bootstrap replications. LL/LP solves,
  resampling laws, inference formulas and already-correct conditional ratios
  are unchanged. Weights that underflow before summation are not recovered.

* Conditional-density and conditional-distribution fallback probes now retain
  the intended physical categorical bandwidths when `bwscaling = TRUE`.
  Valid initial candidates and the existing penalty/restart rules are unchanged.

* External local-linear and local-polynomial fits now report exactly empty
  computed kernel rows as NA in the undefined output components, with one
  warning after the caller finishes. Supported components and required
  training, covariance and bootstrap computations keep their existing policy.
  Point and asymptotic plot calls collect row notices across their evaluation
  grids, and vector-tau least-squares quantile predictions report once per
  completed call. Required bootstrap failures remain terminal.

* Single-index NOMAD refinements retain the original physical bandwidth
  lower bound. A valid incumbent is no longer rejected because its changed
  index coefficients imply a different scale at the refinement handoff.
  Direct fixed-degree calls and actual-beta restart scaling are unchanged.

* Single-index fixed-bandwidth NOMAD searches reject out-of-domain exploratory
  bandwidths without terminating the search. Explicit input and final-result
  validation remain strict, and unexpected callback errors still propagate.

* Generalized-NN single-index training summaries and coefficient covariance
  use a consistent training-evaluation convention when external predictions
  or different inference outputs are requested. Previously affected serial
  external-data results are corrected; the ordinary no-external-data
  convention and MPI counterpart are preserved.

* Shared call evaluation now distinguishes returned condition objects from
  raised errors. MPI worker command loops also accept returned try-error
  values; actual raised-command failure policy is unchanged.

* Coefficient plot controls again validate their values after the dotted-name
  migration. Single-index objects retain their observation count for vector
  indices, and NOMAD-only index searches report cumulative evaluations plus
  the existing final certification count.

* Regression fits and plots replay a saved formula data expression in the
  bandwidth call's existing environment, so wrapper-local data remains usable
  after the wrapper returns and after saving and reloading the bandwidth object.

* Single-index and location-scale quantile NOMAD searches now preserve unexpected
  callback errors instead of returning a surviving candidate. Legitimate typed
  candidate rejection keeps its existing penalty; successful search arithmetic
  and evaluation budgets are unchanged.

* External local-polynomial npreghat mean rows with finite moment systems but
  exactly zero computed kernel weights now return NA with one informative
  warning. Defined rows, kernel accumulation, solver arithmetic and internal
  complete-operator failure behavior are unchanged. The same row contract
  covers matrix, apply, constraint and recomputed prediction outputs.

* Public plot arguments now use dotted names: data.overlay, data.rug,
  factor.boxplot, boxplot.outliers, coef.index, common.scale, proper.method,
  proper.control, boot.control, grid.control and render.control. Plot helpers
  similarly use bar.num, y.vars, y.dat and pair.list. Retired underscore
  spellings are no longer accepted. Function names, saved-object fields,
  defaults and numerical behavior are unchanged; typo hints use the new names.

* Single-index starting-bandwidth restoration now shares the active search
  progress line. Its width-aware notice clears with the search rather than
  leaving a separate message; search computations and failure diagnostics
  are unchanged.

* Finite-data local-constant single-index hat and plot apply helpers now use
  shared kernel moments instead of materializing the full weight matrix.
  Existing bandwidth normalization and denominator floors are preserved;
  explicit matrix outputs and non-finite-input matrix semantics are unchanged.

* Ichimura coefficient covariance now computes its conditional moments and
  common denominator in one kernel-sum call. Covariance formulas, MPI row
  ownership and random-number generation are unchanged; no kernel-weight
  matrix is constructed.

* Single-index fitting now shows one immediate activity line with native
  heartbeats for fitted values and coefficient covariance. Explicit bootstrap
  standard errors retain their separate replication progress, without repeated
  child-regression messages.

* Single-index fits now default to `se = TRUE, se.type = "asymptotic"`,
  providing asymptotic coefficient covariance and fitted-value standard errors;
  `gradients = TRUE` also provides gradients and their asymptotic standard errors.
  Use `se = FALSE` to omit all uncertainty calculations, or explicitly select
  `se.type = "bootstrap"` for the previous fitted/gradient bootstrap standard
  errors. Coefficient covariance remains asymptotic in either mode. Reusing
  `npindex(bws = model$bws, ...)` adds requested outputs without repeating search.
  Singular coefficient-information matrices produce an informative error.
  Local-constant uniform kernels have zero first derivative and cannot provide
  this covariance for free coefficients; use se = FALSE for point estimates.
  Default inference adds training derivative and moment work and can cost
  appreciably more than a point-estimate-only fit when bandwidths are held.

* Adaptive-NN single-index kernel sums now use donor-bandwidth normalization
  consistently for fitted values, bootstrap fits and covariance conditional
  moments, including the shared LC hat/plot helpers. This corrects an existing
  discrepancy between mean-only and gradient fits; fixed and generalized-NN
  weighting is unchanged.

* Single-index pairs bootstraps reuse the already-computed training index,
  avoiding a sampled predictor-matrix copy and repeated matrix-vector product
  per replication. Resampled fits and NN radii are still recomputed. This also
  removes a one-predictor bootstrap dimension-dropping error.

* Single-index bootstrap standard errors now report bootstrap replications,
  forwarding nested regression-fit activity to the same progress display.
  Resampling, numerical results and the number of replications are unchanged.

* Automatic bandwidth calls no longer repeat an entire failed search in other
  calling frames. Original computation errors propagate normally; namespace-only
  calls resolve the package selector before execution.

* Automatic fixed-bandwidth single-index restarts use the index scale of
  their drawn beta, with generated h restricted to the existing search bounds.
  The scale correction can change fresh multistart search results; explicit
  bandwidths, first starts, RNG draws and the search domain are unchanged.
* Kernel-sum formula calls with a subset again allow data to be omitted.
* The internal single-index first-start restoration control is documented.

* Serial smooth-coefficient NOMAD searches restore the first unexpected evaluator
  condition after native completion, before further evaluation or payload/Powell
  work. Typed invalid candidates retain their existing handling.

* Fitted regression releases its owned buffers and partial trees when an error
  leaves the native fit, preserving the original condition and numerical path.

* Explicit na.exclude predictions restore omitted evaluation rows, including
  native evaluation arguments and composed formula/native omission maps.

* LSQ formula dispatch preserves positional subset expressions for the existing formula owner.
* Unconditional NOMAD restart reports and selection now use the existing raw certificate while preserving native recovery scores and official diagnostics.
* Regression MADS restart reports and selection use their existing raw endpoint certificates while preserving native admission and recovery diagnostics.
* Conditional-distribution fixed-degree native restart summaries report their
  existing raw endpoint scores. Official solver diagnostics retain the score
  for the official solution; selection and recovery are unchanged.

* Conditional-density native restart summaries and selection now reuse the
  existing raw endpoint score, retaining native-score admission and recovery
  seed policy. No objective evaluation is added.

* Local-polynomial matrix application keeps training and evaluation bases in
  the same public coordinates when an optional conditioned hat-block attempt
  returns failure. The existing generic computation reuses the retained source
  basis; conditioning, solver and failure policies are unchanged.

* Fixed-bandwidth single-index R optim searches restore raw-invalid automatic
  initial bandwidths by at most eight doublings across all supported kernels,
  orders and regression types, including automatic degree-search children.
  Each restart retains its own index scale. Valid original starts use the same
  optimizer invocation without an extra objective evaluation; rejected first
  scalars are counted and incur no gradient stencil. Explicit or held invalid
  starts and exhausted restoration fail the call. Later trials, NN bandwidths,
  NOMAD's direct search and final certification are unchanged.

* Least-squares quantile formula dispatch preserves named subset expressions for
  evaluation in the formula method's data mask.

* Formula subsets in npcmstest, npqcmstest and npksum resolve data columns
  before caller-local bindings. Data and subset expressions are evaluated
  once.

* Nearest-neighbour bandwidth metadata retains continuous neighbour counts
  with either bandwidth-scaling setting, including semiparametric helpers.
  Conditional bandwidth constructors normalize categorical metadata before
  validation. Fixed continuous bandwidth scaling is unchanged.

* Supplied categorical scale factors are validated in physical kernel units.
  Fixed-bandwidth categorical optimizer transforms now use the corresponding
  scale-factor bounds in both directions; physical kernel caps, penalties
  and optimizer algorithms are unchanged.

* MADS/NOMAD+Powell endpoint checks now evaluate the exact returned
  bandwidths without a round trip through optimizer coordinates. This
  avoids rounding the point being certified; search coordinates, objective
  arithmetic and the number of certification evaluations are unchanged.

* Conditional-distribution NOMAD+Powell refinement compares and publishes
  the raw endpoint objective already obtained by certification, while
  preserving the optimizer's initial-value and history diagnostics.

* NN bandwidth owners release their existing native workspaces before
  reporting constant continuous support. The diagnostic and successful
  search calculations are unchanged; this does not change fitted-value
  error unwinding.

* Regression gradients and their standard errors retain a one-row matrix
  for a single evaluation point, including categorical asymptotic plot
  slices; values and multirow output are unchanged.

* Conditional-distribution cross-validation now rejects failed nearest-neighbour
  bandwidth preparation before using its output. Zero-radius trial candidates
  no longer depend on uninitialized memory; valid-candidate objectives and
  the existing invalid-candidate penalty are unchanged.

* Conditional categorical bandwidth admission now applies the response and
  regressor unordered-kernel bounds independently. Mixed Li-Racine and
  Aitchison-Aitken conditional-density searches no longer reject valid
  response bandwidths or admit response bandwidths above their kernel bound;
  objective arithmetic and same-kernel behavior are unchanged. Selected
  bandwidths and estimates can change when the response and regressor use
  different unordered kernels.

* The embedded wage1 example bandwidth objects were regenerated with current
  constructors. Their printed bandwidths and derived example results can differ
  from earlier releases.

* Cell and exhaustive degree searches now skip only typed LP-admissibility
  and NN-candidate rejections. Unexpected evaluator errors propagate with
  their original condition instead of silently selecting a surviving degree.
  LP rank/capacity rejection messages and successful search arithmetic are
  unchanged; failed searches also release their progress display.

* Conditional bandwidth transforms and categorical baseline probes now
  preserve the raw Y-unordered, Y-ordered, X-unordered, X-ordered layout and
  the corresponding per-side kernel caps. This corrects mistyped coordinates
  when ordered responses and unordered regressors coexist; objective kernels
  and the penalty formula are unchanged.

* Undefined external smooth-coefficient fits now return NA for the affected
  fitted values, coefficients and requested standard errors, with one
  informative warning, instead of publishing an internal finite penalty.
  Valid rows and the ridge schedule are unchanged. Required training and
  bootstrap fits still fail when their local systems cannot be solved.

* Earlier regression uncertainty documentation clarified the then-current
  fitted-residual HC0 calculation, derivatives, contrasts and beta kernels.
  The residual-normalized correction described above supersedes that
  uncorrected formula.

* Formula estimators and regression inference helpers now retain explicitly
  supplied data values in wrapper calls, instead of looking up their argument
  names again in the formula environment, including least-squares quantile
  regression bandwidth and fitted-object formula calls. Formula expressions and retained
  data behavior are unchanged.

* Smooth-coefficient NOMAD bandwidth objects now publish complete bandwidth
  and scale metadata. This fixes multi-smoother summaries and NA fits or
  predictions from malformed scalar degree-zero fixed-bandwidth objects,
  without changing the selected bandwidth, objective or search evaluations.

* Smooth-coefficient bandwidth selection skips redundant scalar zero-ridge
  correction work, retaining the existing ridge safeguards and results.

* Base plots now annotate supported variability intervals consistently for
  levels and gradients/effects: individual full-frame figures retain subtitles
  even with `options(plot.par.mfrow=FALSE)`, while package-managed multi-panel
  pages receive one shared top annotation. Categorical bars and multi-quantile
  overlays use appropriate interval and coverage-scope wording. Explicit
  subtitles and existing graphics geometry are preserved; labels that cannot
  fit safely are omitted. Class-probability plots now honor the same global
  layout option.

* `gradient.order` is the sole supported derivative-order spelling. The
  `gradient_order` plot/single-index alias is no longer accepted, and gradient
  accessors reject that spelling instead of ignoring it. Diagnostics suggest
  `gradient.order` only on routes that support derivative-order requests.

* Mixed fixed/free degree searches now assert after NOMAD returns that every
  fixed degree coordinate still equals its requested value.

* NN recovery now leaves an allowed extended-NN incumbent unchanged and
  proceeds to the existing fail-closed result instead of reporting an internal
  ordinary-domain error.

* Source builds on R versions before 4.4 now provide the LAPACK integer type
  used by compact-support regression checks, while newer R versions retain
  their native LAPACK type definition.

* Long-double entropy-bootstrap work arrays now use R's aligned long-double
  allocator, avoiding undefined behavior on platforms that require stricter
  alignment.

* `npksum(..., permutation.operator = "integral")` now initializes and
  independently owns its tree-range workspaces, preventing invalid memory
  access and double-free behavior without changing the integral result.

* Gradient standard errors for a single-level ordered predictor no longer read
  beyond the predictor's category table, including adaptive-NN dense and tree
  routes; its identically zero contrast remains unchanged.

* Automatic cell-based coordinate/exhaustive polynomial-degree searches now
  preserve a configuration error when every degree candidate fails with the same
  condition. Missing fixed `cker*`, `cxker*`, or `cyker*` limits, and missing
  required beta-kernel support specifications, now produce the same actionable
  messages as NOMAD searches; heterogeneous or numerical candidate failures
  retain the generic no-admissible-model error.

* NOMAD degree searches now honor mixed fixed/free degree bounds such as
  `degree.min = c(0, 0)` with `degree.max = c(1, 0)`. Fixed coordinates remain
  fixed through starts, restarts, optional Powell refinement, and result
  metadata; all-free, all-fixed, and exhaustive searches are unchanged.

* Evaluation-only bandwidth calls now return one initialized history entry,
  rather than exposing unused entries when more than one start was requested.

* Bandwidth searches no longer reuse large-bandwidth geometry from an earlier
  dataset. This fixes history-dependent density cross-validation values and
  protects subsequent searches and regression fits from stale geometry.

* Nearest-neighbour bandwidth searches now certify selected objectives and
  optimizer handoffs against the underlying objective. When an automatic
  search ends on an invalid NN plateau, eligible owners can make one
  deterministic same-optimizer restart from a verified ordinary-NN point.
  The recovery schedule is bounded and does not exhaust all feasible NN
  indices; failure to find a valid restart is not a proof of infeasibility.
  Invalid trial points remain penalized during search; mandatory refinement,
  explicit starting values, and already-valid search paths are preserved.

* Smooth-coefficient NN searches now distinguish categorical smoothing
  bandwidths from continuous NN indices. Single-index searches also preserve
  explicit starting-bandwidth validation instead of silently replacing an
  invalid supplied start.

* The proof-of-concept `npregiv()` and `npregivderiv()` estimators now reject
  unsupported nearest-neighbour bandwidths early with a clear message.
  Fixed bandwidths remain the default and are unchanged.

* Formula-created single-index fits now reuse their retained training data for
  prediction and subsequent evaluation, avoiding spurious missing-variable
  errors when response and predictor columns are not separate workspace objects.

* Terminal nearest-neighbor zero-radius errors now share a diagnostic that
  identifies the continuous variable and tied-donor geometry and suggests
  fixed bandwidths (the default). `npsigtest()` also identifies the failing
  stage. Successful fits, positive radii, and optimizer trial rejection are
  unchanged; repeated values alone are not rejected.

* Principal estimator and bandwidth entry points now reject unknown named
  arguments passed through `...` before evaluating their values. Diagnostics
  suggest an unambiguous canonical spelling when one is available, while
  unnamed arguments and established retired-argument diagnostics retain their
  existing behavior.

* Public bootstrap-count arguments for `npcmstest()`, `npqcmstest()`,
  `npdeneqtest()`, `npdeptest()`, `npsdeptest()`, `npsigtest()`, `npindex()`,
  `npsymtest()`, and `npunitest()` are now consistently named `B`, with the
  historical default of 399. The former `boot.num` input is retired; existing
  `$boot.num` result metadata is unchanged.

* `npsigtest()` now displays one immediate, call-owned progress row with
  compact predictor names. Individual tests report predictor-level completion
  and an ETA after the first test finishes; joint and single-predictor tests
  retain bootstrap-level progress. Statistical results and random-number
  streams are unchanged.

* Eligible individual \`npsigtest()\` calls with fixed bandwidths and IID
  residual resampling now evaluate bootstrap statistics in bounded native
  tiles. The streamed path preserves the incumbent bootstrap samples, random
  number stream, statistics, and P-values; unsupported configurations continue
  through the incumbent R loop.

* Regression, conditional-density, and conditional-distribution evaluation
  rows whose required kernel normalization is exactly zero or non-finite now
  fail closed locally. Their estimate, standard error, gradient, and gradient
  standard error components are returned as \`NA\` where requested, and the
  corresponding hat row is \`NA\`, without invalidating defined rows. Every
  finite nonzero normalization remains literal and unperturbed.

* Fitted density and conditional-density log-likelihood aggregates now use the
  literal logarithm of every represented strictly positive estimate, including
  subnormal values. Exact zero and negative estimates retain the established
  computational placeholder; no underflow-recovery owner is introduced.

* Density likelihood cross-validation now takes the true logarithm of every
  represented strictly positive contribution, including positive subnormal
  values and positive values retained only in signed-log form. Exact zero and
  negative contributions retain the established guarded mapping.

* Leave-one-out regression and conditional objective shortcuts now preserve
  every finite nonzero `1 - h_ii` leverage denominator exactly. Near-saturated
  all-large polynomial designs no longer replace a positive representable
  denominator by machine epsilon; zero or non-finite deletion denominators
  cede to the established invalid-objective path.

* Continuous generalized- and adaptive-nearest-neighbour bandwidths now use
  the literal kth admitted distance throughout ordinary, beta-kernel, and
  kernel-sum owners. Exact ties that make this radius zero invalidate the
  requested bandwidth instead of silently widening it to the nearest
  positive distance or replacing it by a numerical floor. Search candidates
  retain their established invalid-objective mapping, while direct public
  computation reports the zero-radius condition. Every represented strictly
  positive radius remains admissible; k, occurrence exclusion, extended-NN
  scaling, cache policy, and public interfaces are unchanged.

* Mixed beta-kernel conditional hats now pass the validated categorical-
  compression state to their private direct kernel-weight owner. This repairs
  conditional density gradients and simultaneous beta-X/beta-Y fits that
  previously stopped for both dense and compressed execution; the two modes
  remain numerically identical.

* Copula density helpers now translate the public normalized ordered
  Li--Racine kernel name to the private `kbandwidth`/`npksum` spelling. This
  restores ordered and mixed copula density semantics without changing public
  kernel names or defaults.

* Wide fixed-bandwidth local-polynomial `npreghat(..., output = "apply")`
  calls now use the training-coordinate block shortcut only when evaluation
  coordinates actually equal the training coordinates. Equal row counts no
  longer misclassify a distinct external evaluation grid; the cutoff,
  arithmetic owners, and public interface are unchanged.

* Partial-linear INID bootstrap coefficient solves now attempt the identified
  unregularized weighted least-squares system before entering their bounded
  ridge ladder. The former positive `1e-12` default altered every identified
  replicate and could make an exactly singular system appear solved before
  the declared `1/n` stabilization step. Failed zero-ridge solves still enter
  the existing `1/n`, `2/n`, ..., capped sequence; public interfaces are
  unchanged.

* Partial-linear bootstrap coefficient solves now ridge the
  Robinson-residualized regressors symmetrically. The former first-coordinate
  RHS adjustment treated a user regressor as an intercept and made ridged
  results depend on regressor order; weights, ridge ladder, resampling, and
  public interfaces are unchanged.

* Smooth-coefficient fit and bandwidth ridge corrections now use the finite,
  nonzero signed pristine moment intercept. This preserves constant
  reproduction under common moment scaling without changing ordinary
  zero-ridge solves or the existing bounded ridge sequence.

* Ridged local-polynomial response, influence, and R fallback solves now
  restore the intercept from the signed pristine Gram intercept rather than
  the already-ridged diagonal or a numerical floor. This restores constant
  reproduction on rank-deficient systems without changing ridge admission,
  sequence, cap, diagnostics, or ordinary zero-ridge arithmetic.

* Regression bandwidth search now preserves terminal invalidity through
  finalization and raw-certifies the selected point once outside the optimizer
  loop. Finite invalid guidance uses the exact criterion-specific constant
  null for CVLS, the trace-one null for CVAIC, and the Bernoulli null for
  CVKS; penalty multipliers must be finite and strictly greater than one.

* Location-scale quantile bandwidth selection now evaluates leave-one-out
  check loss against the requested response when its fitted transformed
  response differs, uses the exact constant-quantile check-loss null for
  invalid guidance, and rejects an optimizer payload whose selected point is
  raw-invalid. The chosen NOMAD/Powell semantics are unchanged.

* Single-index bandwidth selection now raw-certifies selected Ichimura and
  Klein--Spady candidates and uses their exact declared leave-one-out nulls
  for outer invalid guidance. Inner regression evaluation returns raw
  objective values or terminal invalidity and propagates implementation
  errors; valid objectives, caches, and optimizer choices are unchanged.

* Positive-degree local-polynomial categorical gradients/effects now report
  their documented fitted-value endpoint contrasts across regression,
  conditional density/distribution, quantile regression, conditional mode,
  location-scale quantile regression, and plot/bootstrap workflows instead of
  native zero placeholders. Categorical asymptotic effect standard errors
  remain unavailable and are reported as `NA`.

* `npsigtest()` now chooses pivotal continuous and unstandardized categorical
  statistics by default, records the effective choice, and fails explicitly
  when requested pivot standard errors do not define a statistic. Bootstrap
  values equal to the observed statistic are included in the empirical upper
  tail, so a completely smoothed-out predictor returns `P = 1` rather than a
  false rejection. Explicit local-polynomial tests require degree at least one
  only for tested continuous coordinates. Joint/predictor status now stays on
  the existing live progress line instead of leaving permanent console notes.

* Multi-response `npreghat(..., output = "apply")` results now retain the
  response column names consistently across native, leave-one-out, and
  row-local local-polynomial owners, matching the corresponding hat-matrix
  application without changing numerical values or ridge diagnostics.

* Local-polynomial response, influence-row, regression-CV, conditional X-row,
  and bootstrap owners now share one typed numerical-rank and ridge policy.
  Ordinary full-rank response rows retain the one-call `DGESV` lifecycle;
  only an ambiguous working-precision pivot enters a cold, symmetrically
  equilibrated singular-value check. Rank-deficient rows use one
  Gram-relative bounded-ridge transcript for response and adjoint solves.
  This restores dense/tree objective parity for sparse fixed-bandwidth LP
  windows and makes `npreghat()` matrix/apply results agree with fitted values
  on deficient rows without adding evaluation-by-training storage.

* Public `regtype="ll"` is now exclusively an API and reporting alias for the
  canonical raw GLP degree-one engine. Regression fit, uncertainty, objective,
  `npreghat()` matrix/apply/constraint, single-index hat, conditional X-side,
  and plot/bootstrap consumers all select owners from canonical engine
  metadata. The retired private local-linear matrix/ridge and single-index hat
  routes have been removed. Results can change only where those private routes
  formerly selected a different regularized fit; explicit raw LP degree one
  and unaffected LC/LP configurations retain their canonical results.

* Fixed local-polynomial plot bootstraps now solve count-compressed response
  moments through the same canonical bounded-ridge response owner used by
  ordinary regression fits. Regression, single-index, smooth-coefficient,
  partial-linear, and conditional density/distribution bootstrap consumers
  share one native batched response solve; the former private determinant,
  `solve()`, and absolute-ridge implementations have been removed. Count
  construction, moment arithmetic, chunking, progress, and MPI orchestration
  are unchanged. Results can change for systems where the retired private
  path imposed a different ridge decision.

* Eligible fixed and generalized-nearest-neighbour LP mean applications now
  use the canonical regression response owner for one or more right-hand
  sides. A one-column `npreghat(..., output = "apply")` no longer falls
  through a separate kernel-weight matrix reconstruction that could disagree
  with `npreghat(..., output = "matrix") %*% y` and `npreg()` on higher-degree
  or low-rank local windows. Adaptive and extended nearest-neighbour routes,
  beta kernels, derivative operators, ridge policy, and public interfaces are
  unchanged.

* Ordinary univariate generalized-nearest-neighbour regression now admits
  `k = 1` only for the mathematically certified positive scalar capability:
  second-order Gaussian or Epanechnikov kernels with LC mean/derivative
  operation, including the equivalent generalized-LP degree-zero
  representation. One exact-field capability resolver supplies manual
  validation, fixed-degree optimizer bounds, and the native lower endpoint;
  an occurrence-aware geometry validator rejects literal zero radii. Adaptive
  NN, mixed or multivariate regression, positive-degree LP (including the
  public `ll` alias), other kernels/bounds, automatic degree domains that can
  reach positive degree, and semiparametric owners retain `k >= 2`. The policy
  contains no sample-size, timing, host, rank-count, or acceleration branch and
  introduces only linear-in-evaluation geometry scratch.

* Regression leave-one-out hat-matrix and apply routes now share the objective's
  occurrence-aware generalized/adaptive NN geometry, source basis, weighted
  moment layout, factorization, bounded ridge transcript, and intercept
  restoration. Matrix and apply consume the same completed influence row;
  apply retains linear-in-sample auxiliary storage and does not construct an
  evaluation-by-training matrix. The reconstructed owner preserves the
  August-16 one-call response solve and existing bounded ridge boundary;
  all-zero compact-kernel LOO rows now fail explicitly instead of returning a
  fabricated ridged estimate.

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
  invocation-owned outer-product workspace across full-contiguous tree rows.
  The unchanged design matrix is packed once while response/weight scratch is
  reused; sparse, partial, and multi-segment support retains the incumbent
  row-local owner. This removes repeated large workspace allocation and
  design packing without changing public bandwidths, objectives, evaluation
  accounting, or fitted results.

* `nplsqregbw()` and `nplsqreg()` now support the package-wide
  `nomad = "auto"` shortcut. One continuous smoothing dimension uses the
  exhaustive degree lattice and jointly optimizes bandwidth and `delta` at
  every degree; two or more dimensions use the existing NOMAD/Powell joint
  search. Historical `nomad = FALSE` and `nomad = TRUE` routes are unchanged,
  and the separate `nomad.pilot` control remains logical-only.

* Conditional-density automatic degree search now provisions response-side
  categorical workspace for every engine that a later degree candidate can
  activate, rather than only the engine selected by the initial degree. This
  restores `npcdensbw()` and `npconmode()` search-trajectory parity across
  serial and MPI owners without changing fixed-candidate objective arithmetic,
  degree policy, or public evaluation accounting.

* On Apple silicon, fixed-bandwidth compact-support regression-tree CV
  objectives of local-polynomial basis width six now use a compile-time
  resident sparse-pair microkernel. It preserves tree and pair order, the
  unique upper Gram triangle, and the shared mirror/solve finisher while
  vectorizing only the contiguous moving-row update. Widths below and above
  six are unchanged, and portable builds retain the incumbent generic width-
  six owner.

* Fixed-bandwidth compact-support regression CV candidates whose kernel
  support covers the complete observed continuous range now use the canonical
  dense LP owner instead of unproductive full-tree traversal. A shared,
  operator-aware geometry classifier also routes full-support fixed kernel
  sums to their dense sibling and searches mixed-support trees only over
  coordinates whose kernels can prune, while continuing to evaluate every
  kernel factor. These exact, allocation-free policies apply to Apple and
  portable builds without sample-size or runtime thresholds; width-two
  regression CV, genuinely sparse, nearest-neighbour, non-compact, and
  unsupported owners retain their established routes.

* Conditional-density native degree search now constructs its finite invalid-
  candidate penalty once from the same canonical baseline-and-retry rule as
  the ordinary bandwidth route, retains that value in the prepared context,
  and restores it before every candidate. This removes candidate-dependent
  penalty recomputation, preserves public evaluation accounting, and carries
  degree-search parameters through the retry workspace without quadratic
  sample storage.

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

* Fixed-bandwidth degree-zero conditional-density CVML with two or three
  compact-support explanatory dimensions now retains sparse tree traversal in
  the canonical scalar joint/marginal owner instead of entering the dense LP
  block stream. The selector is based only on prepared engine topology and
  preserves `np.tree = "auto"` kernel eligibility, LC/raw-LP0/Bernstein-LP0
  equivalence, linear auxiliary storage, and existing wider, positive-degree,
  nearest-neighbour, dense, and MPI owners.

* `npconmode()` now evaluates categorical response levels through one
  memory-bounded conditional-density batch instead of rebuilding the same fit
  lifecycle for every level. Requested class-probability gradients retain one
  singleton evaluation for their selected level while all other levels share
  the batch. The canonical `npcdens()` engine, probability repair, class
  selection, and result arithmetic are unchanged; temporary evaluation rows
  are capped independently of the number of response levels.

* Eligible fixed-bandwidth local-polynomial CV objectives now apply the same
  basis-neutral owner rule to generalized, additive, and tensor bases:
  resident accumulation through basis width six, packed BLAS from width seven,
  and the independently qualified generic resident owner for sparse tree rows.
  The hoisted basis, canonical solve, deletion, projection, and loss arithmetic
  are shared; no basis-family or kernel-spelling fallback remains. Widths five
  and six retain only the unique symmetric Gram triangle during accumulation
  and mirror it once before solving.

* Fixed-bandwidth compact-support `npscoef` objectives now pack each sparse
  tree row's active support into the existing BLAS outer-product engine for
  moderate and wide local-polynomial moment systems. Dense, degree-zero,
  nearest-neighbour, signed, derivative, and permutation routes retain their
  existing arithmetic and storage contracts; scratch is bounded by
  `support * (width_A + width_B) + width_A * width_B`, never by sample pairs.

* Fixed-bandwidth degree-zero conditional-density CVLS now recognizes an
  all-categorical explanatory side at its constant upper kernel limit as the
  canonical global-X design. This restores the unconditional-density limit,
  including the ordered Li-Racine case, and avoids redundant conditional X
  weighting. Generalized/adaptive NN and response-side kernel limits are
  unchanged.

* Degree-zero conditional-density and conditional-distribution CVLS objectives
  now enter the same canonical all-large explanatory-bandwidth context as
  positive-degree local polynomials.  The width-one case uses its implicit
  unit basis and scalar Gram inverse, avoiding the dense locally weighted pass
  while preserving delete-one objective semantics and linear auxiliary
  storage. Conditional-distribution tree rows now translate their independent
  response permutation once at the shared indicator boundary.

* Fixed-bandwidth two-predictor local-linear CV objectives now keep the sparse
  tree row resident and accumulate the same unique symmetric Gram triangle as
  the canonical dense width-three engine before using their shared mirror,
  solve, and uncentered projection. Compact-support tree support and kernel
  semantics are unchanged; the redundant lower-triangle pair updates are
  removed without additional sample-sized storage.

* Bootstrap-intensive tests now reuse fixed kernel geometry in bounded
  replication groups. `npcmstest()` and `npqcmstest()` keep one model refit
  per draw but contract the resulting residual or score columns through
  multi-response kernel sums; fixed, unbounded `npdeneqtest()` routes process
  pooled multiplicity columns while bounded and nearest-neighbour routes keep
  literal duplicate resampling; and `npdeptest(method = "summation")` uses a
  bounded native index batch that shares only the invariant response-side
  marginal work. Draw order, complete bootstrap payloads, bandwidth and model
  semantics, RNG restoration, progress, and linear auxiliary-memory scaling
  are preserved.

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
  `npsdeptest()`. Univariate and symmetry bootstraps consume bounded chunks of
  multiplicities while preserving the exact incumbent `sample.int()` or
  `boot::tsboot()` draw plan and RNG state; dependence statistics fuse their
  marginal and joint kernel sums without an n-by-n cache. Explicit
  generalized- and adaptive-nearest-neighbour requests remain on literal,
  duplicate-preserving resample/refit routes so sample-owned radii and ordered
  lag pairs are unchanged.

* General bounded conditional-density CVLS now reuses each response-side
  quadrature kernel tile across a memory-bounded group of explanatory-row
  tiles. Ordinary and beta kernels share the same canonical engine; objective
  arithmetic and accumulation order are unchanged, and storage remains
  `O(n B)` under the existing conditional-LP tile budget.

* Eligible bivariate surface plots now share one base-perspective frame owner,
  one viridis palette mapper, and one neutral-gray three-face grid style across
  the base and `rgl` renderers. The shared renderer is used by the regression,
  density, distribution, conditional, semiparametric, conditional-mode, and
  copula plot families; estimator and plot-data arithmetic are unchanged.

* Uncertainty migration: `npudens()`, `npudist()`, `npcdens()`,
  `npcdist()`, `npplreg()`, `npqreg()`, `npconmode()`, and `npcopula()`
  now default to `se = FALSE`, joining the existing opt-in controls in
  `npreg()`, `npscoef()`, and `nplsqreg()`. This is intentionally
  backward-incompatible for scripts that assumed these uncertainty fields
  were always populated in earlier releases or candidate builds. The
  exception is `npindex()`: its default `se = TRUE` includes asymptotic
  coefficient covariance; explicit `se = FALSE` omits it. Gradients and
  uncertainty remain separate requests where supported. For example,
  `npreg(bws = model$bws, gradients = TRUE, se = TRUE)` adds requested
  outputs without repeating bandwidth selection; supply the original training
  data if it cannot be recovered from the bandwidth object. Here `model`
  is an illustrative object name. `se()`, `vcov()`, `coef(..., se = TRUE)`
  and `gradients(..., se = TRUE)`, where supported, only extract stored
  results and give a no-search refit message when those results are absent.
  Explicit inference requests such as `predict(..., se.fit = TRUE)` and
  `plot(..., errors = "asymptotic")` request their needed computation.
  Omitted internal uncertainty fields can be `NULL` or length-zero vectors
  (for example, `derr = numeric(0)`); use the public extractors rather than
  relying on one internal empty-field shape. The former Boolean estimation
  argument `errors` remains rejected; plot methods retain their method-valued
  `errors` control. No bandwidth search is repeated merely to extract stored
  results.

* Conditional density and distribution kernel summaries now identify
  explanatory and dependent kernels separately for continuous, unordered, and
  ordered variables. Bandwidth-object `print()` and `summary()` methods share
  one formatter, public output no longer exposes internal `cxker*`/`cyker*`
  argument mnemonics, and invalid categorical-kernel selections name the
  specific `uxkertype`, `uykertype`, `oxkertype`, or `oykertype` argument.

* Interactive `renderer = "rgl"` surfaces now use a high-DPI-aware widget
  canvas, a higher-resolution uncertainty legend, and an initial camera more
  closely matched to the base perspective. The renderer caps backing-store
  scaling at two device pixels per CSS pixel, retains native rgl axis labels
  and interaction, and preserves all existing `rgl.*` user overrides. Shared
  rgl setup and option routing now have one internal owner across supported
  plot families. Existing `band = "all"` legends now use the same
  construction-aware quantile, rank, and normal-SE labels as the base renderer,
  without adding an rgl subtitle or footer; estimator, bootstrap, interval,
  and base-renderer arithmetic are unchanged.

* `npindexbw()` now canonicalizes its numeric design matrix once at the public
  boundary and reuses that matrix during internal index and coordinate setup.
  This removes redundant full-matrix copies before objective evaluation while
  preserving estimator arithmetic and conditions.

* The integration routes for `npunitest()`, `npsymtest()`, `npdeptest()`, and
  `npsdeptest()` now use bounded deterministic corrected-trapezoid quadrature
  for their default Gaussian statistics. Univariate iid bootstraps reuse a
  fixed grid and evaluate bounded chunks of resampling counts; symmetry
  bootstraps reuse reflected density values; and bivariate quadrature calls
  the registered constant-storage Gaussian evaluator in bounded tiles. This
  removes the `cubature` dependency while preserving bootstrap sampling laws,
  fixed bandwidth reuse, result shapes, seed restoration, and nondefault
  kernel routes.

* Removed the private, unexported package-local `gsl.bs` bridge and its
  registered native wrappers. Public spline construction remains owned by
  `crs::gsl.bs()`; the independent B-spline primitives used by the canonical
  Bernstein local-polynomial engine remain in place, with estimator arithmetic
  unchanged.

* The canonical beta row engine now owns its malloc-backed route, categorical,
  derivative, and row workspaces behind one unwind-protected cleanup boundary.
  Nested derivative helpers borrow that invocation-scoped scratch rather than
  allocating it per call, so user interrupts and R errors cannot strand native
  storage. Successful PDF, CDF, convolution, derivative, powered-derivative,
  and regression arithmetic is unchanged.

* Beta nearest-neighbour bandwidth preparation now calls the shared
  continuous-distance owner directly. It no longer temporarily rewrites
  unrelated package-global scaling and categorical factors, while retaining
  the same generalized- and adaptive-nearest-neighbour arithmetic.

* Nearest-neighbour lookup and distance failures now leave both shared
  bandwidth owners through their normal cleanup exit. This removes retained
  standard-deviation and distance workspaces after recoverable failures while
  preserving successful fixed, generalized-nearest-neighbour, and adaptive-
  nearest-neighbour arithmetic and conditions.

* Removed the unreachable conditional-bootstrap exact-state sidecar and its
  dead-state-only test assertions. Public conditional bootstrap routes retain
  their existing canonical refit, count, and fixed local-constant owners.

* Removed the unregistered and callerless internal `beta_kernelsum` facade.
  Beta bandwidth-mode ownership now resides in the shared bandwidth contract;
  all estimator, kernel-sum, LP, and optional-output routes continue to use the
  canonical continuous-row engine.

* Higher-order beta convolution now prepares each distinct side/component
  shape once per observation pair rather than rebuilding it for every
  opposing component. The stack-local decomposition retains the canonical
  centered gamma-ratio and signed component-pair accumulation order, changes
  no order-two arithmetic, and adds no dynamic or sample-sized storage.

* Fixed-bandwidth beta CDF rows now prepare each coordinate/component
  concentration once per invocation rather than recomputing it for every
  evaluation-observation pair. The fixed-only sibling shares the canonical
  CDF arithmetic and failure statuses with scalar and nearest-neighbour rows,
  cannot fall back after selection, and adds only coordinate-by-component
  transient storage.

* Beta CDF rows now prepare observation support coordinates and higher-order
  coefficient state once per invocation rather than rebuilding them for every
  evaluation-observation pair. Fixed, generalized-nearest-neighbour, and
  adaptive-nearest-neighbour distribution rows share the same canonical
  `pbeta()` arithmetic and failure contract. The prepared route is selected
  outside row traversal, cannot fall back to the scalar CDF evaluator, and
  uses only beta-coordinate-linear transient storage; density-only routes
  retain their existing prepared-PDF path unchanged.

* Beta regression gradient rows now form the signed-log level and its
  regular/jump target derivative from one canonical component preparation.
  The level phase still completes before derivative work, preserving native
  failure precedence and exact fitted, gradient, and standard-error
  arithmetic while removing the duplicate beta shape and PDF pass. Fixed- and
  generalized-nearest-neighbour gradient rows additionally reuse their
  prepared evaluation-row components. Adaptive-nearest-neighbour gradient rows
  reuse invocation-owned observation transforms and coefficient state while
  retaining pair-owned bandwidth shapes and normalizers. Each bandwidth
  topology selects an isolated sibling before row traversal, with no inner-loop
  route branch or scalar fallback. Derivative-only fixed/GNN state is allocated
  lazily, so level-only fits retain their incumbent preparation cost.

* Adaptive-nearest-neighbour beta PDF rows now prepare support-transformed
  observations and higher-order coefficient state once per invocation while
  retaining observation-owned bandwidth shapes and normalizers in the
  pairwise evaluator. This preserves exact kernel-weight, objective, fit, and
  prediction arithmetic, uses beta-coordinate-linear transient storage, and
  cannot fall back to the scalar row after the prepared route is selected.

* Fixed- and generalized-nearest-neighbour beta PDF rows now prepare support-
  transformed observations, component normalizers, and higher-order
  coefficient logs at their invocation or evaluation-row lifetime instead of
  recomputing them for every observation pair. The shared row engine retains
  exact objective, fit, prediction, and kernel-weight arithmetic, has no
  runtime fallback after route selection, and uses only beta-coordinate-linear
  transient storage.

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
  the scalar width-one or general LP engine through the shared canonical rule.
  Returned objectives continue to replay exactly without fallback or
  post-search rewriting.

* Derivative kernel-weight output now uses invocation-scoped ownership and
  each derivative block's canonical tree-support intervals. This removes
  mutable native output state and corrects compact-support adaptive-nearest-
  neighbour tree output while preserving exact dense, tree, automatic, and
  repeated `npksum()` results, including mixed score/OCG blocks.

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

* Corrected conditional plot/bootstrap dispatch so nonfixed local-polynomial
  beta fits retain the declared LP engine on every exact resample instead of
  being intercepted by the local-constant count evaluator.

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
  objective from mixed legacy state.

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
  checked 64 MiB workspace. Above the ceiling, stable resident rows and
  dedicated scratch rows are filled through the shared categorical-profile
  kernel engine. A 32-profile convolution supertile preserves every scalar
  accumulation order while reusing each response-convolution row across the
  group, reducing representative objective time without changing objective
  bytes. Internal failure after route commitment is terminal rather than a
  silent fallback.

* Fixed-bandwidth categorical conditional-distribution CVLS now bounds its
  aggregate predictor- and response-profile kernel workspace at 64 MiB.
  Problems that fit retain the established dense arithmetic exactly; larger
  profile topologies use stable reuse-priority resident rows plus one scratch
  row, filled through the shared categorical-profile kernel engine. Objective
  accumulation order and values are unchanged, and a failure after route
  commitment can no longer silently enter a different implementation.

* Fixed-bandwidth ordered-categorical unconditional-distribution CVLS now
  consumes compressed kernel rows through a checked, caller-owned tile capped
  at 64 MiB instead of retaining the full evaluation-by-training profile
  matrix. The traversal and loss accumulation orders are unchanged, and
  objective values remain byte-identical; the general route remains canonical
  when a single row cannot fit the bounded workspace. Once the bounded route
  is selected, internal failure now fails cleanly rather than silently
  re-entering the general route.

* Added a dormant, rank-local categorical-profile tile engine for subsequent
  migration of dense profile consumers. Its caller-owned output is bounded by
  a checked 64 MiB ceiling, immutable training profiles can be validated once
  per traversal, and no production estimator dispatch changes in this
  checkpoint. The same work corrects unordered derivative dispatch so score
  operators reach their canonical kernels instead of the historical normal-
  kernel fallback; ordinary objective, fit, and prediction results are
  unchanged.

* Powell-side objective caches now use checked capacity, key-width,
  load-factor, and rehash arithmetic under an independent 64 MiB peak ceiling
  per native table. If a table cannot grow, its existing entries remain
  available while the selected optimizer continues without further insertion.
  Cache keys, objective values, optimizer behavior, and the strict
  `np.objective.cache` option contract are unchanged.

* Smooth-coefficient local-polynomial cross-validation and fitting now solve
  the common stable zero-ridge row systems through one registered native
  entry, reusing a bounded `O(p^2 + p)` LAPACK workspace instead of making one
  R-to-LAPACK transition per evaluation row. Successful wider batches also
  project their hoisted basis rows in ascending basis order through one
  register-local native entry instead of slicing, reshaping, and dispatching
  every row in R. Width-one and width-two explicit solvers are unchanged; if
  any wider row is non-finite, singular, or ill-conditioned, the entire batch
  returns to the established R ridge loop. Estimator formulas and ridge policy
  are unchanged, although LAPACK workspace alignment and projection
  reassociation can change last floating-point bits.

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
  response-plus-basis operand once per objective evaluation for eligible
  full drop-one rows, instead of repacking it for every observation.
  Adaptive-nearest-neighbour, tree, reduced-row, and specialized resident-row
  routes are unchanged, as are objective values.

* Conditional-density and conditional-distribution local-polynomial
  cross-validation now retain signed higher-order predictor-kernel weights in
  every delete-one X row. The obsolete QR route silently discarded negative
  weights; it has been removed in favor of the canonical signed full-row solve
  and exact diagonal deletion, with a sign-preserving denominator. This
  intentionally corrects affected higher-order CVLS and CVML objectives;
  ordinary second-order results are unchanged.

* Conditional local-polynomial full-row deletion now preserves every finite,
  nonzero signed `1 - H_ii` denominator exactly rather than replacing
  sub-machine-epsilon values by a fixed floor. Exactly zero or non-finite
  denominators fail through the existing clean objective path because the
  corresponding deleted system has no valid finite row. Ordinary
  non-singular objectives are unchanged.

* Eligible fixed and generalized-nearest-neighbour conditional
  local-polynomial full-row blocks now reuse the bounded weighted-design
  BLAS assembly already used by conditional-density CVLS. This accelerates
  conditional-density likelihood cross-validation and conditional-
  distribution least-squares cross-validation without changing the scalar
  width-one, adaptive-nearest-neighbour, non-Accelerate, or allocation-
  fallback paths. The formulas are unchanged but BLAS reassociation can
  change last floating-point bits.

* Adaptive-nearest-neighbour conditional local-polynomial objectives now
  reuse the same bounded weighted-design BLAS assembly when the basis has at
  least four terms. Widths one through three retain their established scalar
  transcript, as do small-sample, non-Accelerate, and allocation-fallback
  routes. This accelerates conditional-density likelihood and least-squares
  cross-validation and conditional-distribution least-squares
  cross-validation; the formulas are unchanged but BLAS reassociation can
  change last floating-point bits.

* Eligible adaptive-nearest-neighbour conditional-density objectives now
  cache ordinary-Gaussian response-bandwidth reciprocals in their existing
  row context, replacing repeated vector divisions with multiplications.
  Non-Gaussian, higher-order, bounded, mixed-response, tree, non-Accelerate,
  and allocation-fallback routes retain the established implementation. The
  formulas are unchanged but reciprocal multiplication can change last
  floating-point bits.

* On Apple ARM64, eligible fixed and generalized-nearest-neighbour products of
  two or more ordinary Gaussian continuous kernels now use one vector
  exponential of the summed squared standardized distances. Non-Apple,
  one-dimensional, tree, bounded, higher-order, adaptive-nearest-neighbour,
  convolution, score, permutation, and generalized-NN large-bandwidth
  shortcut routes retain the established implementation. The formulas are
  mathematically identical but can differ in their last floating-point bits.

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
  systems fail instead of retrying indefinitely, and the fit route releases
  its native workspace before reporting failure; ordinary successful systems
  retain their established solve and ridge sequence.

* Aligned the dedicated `plot()` methods for `npregiv` and
  `npregivderiv` with the package-wide regression plotting vocabulary. They
  now use `gradients`, `data_overlay`, and `data_rug`; `npregiv` defaults to
  its structural level with the training response overlaid, while
  `npregivderiv` continues to default to its derivative. Automatic limits
  include every active curve, response overlay, and training-support rug,
  without replacing available fit-time evaluation curves. The experimental
  plot-only controls `plot.data`, `deriv`, and `phi` have been removed and now
  fail with migration guidance. Estimator objects and numerical results are
  unchanged.

* Conditional density and distribution level plots now support bootstrap
  intervals when the continuous X and Y sides use different kernel families
  or orders, including beta on either side. Beta-containing bootstrap levels
  use the same signed log-domain numerator/denominator accumulation as the
  fitted estimator; this also corrects same-family higher-order beta
  bootstrap values when a finite negative or extremely small explanatory
  kernel sum was previously replaced by a positive machine-epsilon guard.
  Fixed, frozen generalized/adaptive nearest-neighbour, and exact
  nearest-neighbour refit semantics are covered. Matching legacy-kernel paths,
  estimator fits, gradients, bandwidth selection, and public defaults are
  unchanged.

* Aligned `npregiv()` and `npregivderiv()` summary output with `npreg()`.
  Fixed smoothing now uses the canonical `Kernel Regression Estimator:`
  description (`Local-Constant`, `Local-Linear`, or the local-polynomial
  degree and basis details) instead of separately printing the legacy `p`,
  regression type, and degree fields. This is a reporting-only change;
  fitted objects, smoothing choices, and numerical results are unchanged.

* `npksum()` now validates its public logical controls and `kernel.pow` before
  native dispatch. This prevents malformed scalar controls from being masked
  by short-circuit evaluation and prevents empty, missing, non-finite, or
  vector-valued kernel powers from reaching C; an empty kernel power could
  previously terminate the R process. Documented scalar inputs and valid
  integer powers retain their established behavior.

* Repaired `npksum()` score and OCG output for mixed data. `compute.score=TRUE`
  now allocates its native result buffer, and packed score, OCG, and continuous
  permutation blocks (including derivative kernel weights) are returned in
  original data-column order. Categorical base-kernel factors are now retained
  in mixed-data continuous permutation blocks, and score-only ordered kernels
  no longer dereference OCG-only state. The defects could previously cause a
  process-level crash, fail during R reconstruction, or return an incomplete
  mixed-data derivative.

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
  bandwidth divisors. Univariate Tikhonov results are unchanged.

* Separated `npregivderiv()` training and evaluation state. A fit-time
  `zeval`/formula `newdata` grid now affects only `phi.prime.eval` and its path;
  training operators, residuals, centering, stopping, selected state, fitted
  values, and training gradients remain unchanged. Arbitrary positive grid
  sizes are supported, while a different `weval` is rejected because the
  inverse problem is defined on the training instruments. The default Issue
  57 local-linear training route is numerically unchanged.

* Hardened the IV public surface: ordinary IV regression stages now own
  `bandwidth.divide`, `ukertype`, and `okertype` without duplicate-argument
  failures and use kernel-appropriate categorical search bounds;
  `npregivderiv(random.seed=)` now controls its internal bandwidth searches;
  scalar controls fail early; summaries report continuous and categorical
  counts and derivative residual bandwidths; and `fitted()`, `residuals()`,
  `gradients()`, and plot methods consistently distinguish training from
  evaluation fields. Post-fit `predict()` remains deliberately unavailable.

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
  stopping behavior, and smoothing choices are unchanged.

* Modernized the public `npregiv()` and `npregivderiv()` interfaces while
  preserving their established numerical defaults. Both now support explicit
  IV formulas (`y ~ z | w`, with optional `| x` where the estimator supports
  it), `data`, `subset`, `na.action`, fit-time `newdata`, fixed
  `regtype = "lc"`/`"ll"`/`"lp"` and scalar `degree` controls, structured
  summaries, and `fitted()`, `gradients()`, and training-row `residuals()`.
  Objects retain compact bandwidth, smoothing, and stage metadata. Legacy
  `p` remains supported by `npregiv()` with unchanged default behavior.
  Unsupported post-fit prediction, derivative `x`, and automatic-degree
  NOMAD routes are documented and fail explicitly instead of being silently
  accepted or approximated.

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
  residual. This can change derivative trajectories, stopping states, and
  fitted curves. The `npregiv()` and `npregivderiv()` examples now both use
  `n = 500`, for which the documented seed produces stable estimates.

* Corrected Landweber-Fridman state indexing in `npregivderiv()`. Iteration
  `N` now consistently denotes `N` completed derivative updates: column `N`
  of `phi.prime.mat` and `phi.mat`, `norm.stop[N]`, `num.iterations`, and the
  returned derivative/function now identify the same state. The initialization
  remains state zero, and evaluated stopping-rule overshoots remain available
  in the iteration matrices while `num.iterations` identifies the selected
  state. This can change the selected iteration and estimates because the
  stopping criterion is now evaluated as `N` times the residual norm at state
  `N`, as described by Florens, Centorrino, and Racine.

* Corrected the empirical adjoint in `npregivderiv()` to use the ordinary
  kernel CDF required by Equation (14) of Florens, Centorrino, and Racine.
  The integral kernel sum had retained an extra continuous-bandwidth factor,
  which could make the Landweber-Fridman stopping norm increase and bias the
  recovered conditional mean. The private adjoint now owns its required
  normalization, while public `npksum()` defaults and regression argument
  forwarding remain unchanged. The iteration guard also examines only the
  computed prefix of the preallocated stopping vector. This resolves issue
  #57.

* Retired the unused experimental truncated-Gaussian continuous kernel and
  its `nptgauss()` configuration helper. The supported continuous kernels are
  Gaussian, Epanechnikov, and uniform; their public interfaces, native codes,
  and numerical behavior are unchanged.

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
  derivatives, fitted values, standard errors, and the native evaluator use
  the same term definition. Local constant, local linear, univariate,
  additive, tensor, and non-local-polynomial semantics are unchanged.

* Hardened native NOMAD observer and interrupt handling. Explicit user
  interruption is now reported as an R `interrupt` condition only after
  native cleanup, ordinary observer errors remain fail-open, native callables
  are resolved per solve rather than retained across package reloads, and the
  declared R 3.5 compatibility floor is preserved.

* Restored timely progress updates during native NOMAD bandwidth and degree
  searches, including long compiled objective evaluations. Iteration, current
  degree, and accepted-best details now follow the package-wide progress
  interval without changing optimization results or evaluation accounting.
  This requires `crs` 0.15-46 or later.

* Restored partial mixed-degree local-polynomial gradient evaluation and
  plotting across regression, conditional density/distribution, least-squares
  quantile regression, conditional quantiles, and conditional modes. A
  continuous derivative requested above that coordinate's fitted degree is
  retained as `NA` and shown as an empty panel while other supported
  derivatives and categorical first-difference effects continue to plot.
  Requests with no available component now fail early with the affected
  predictor orders and degrees and actionable refitting guidance.

* Corrected common-scale uncertainty-band rendering for partially available
  gradient plots. Empty derivative panels now retain their list positions and
  finite-range calculations no longer emit spurious `min`/`max` warnings;
  supported estimates and interval bands are unchanged.

# np 0.70-5

* Corrected the heteroskedasticity-robust Ichimura index-coefficient
  covariance returned by `vcov(npindex(..., gradients = TRUE))` when the model
  has more than one free index coefficient. Link-gradient and residual weights
  are now applied observation by observation. Fits with one free coefficient
  remain exactly unchanged, as do beta, bandwidth, objective, fitted-value,
  residual, and gradient results.

* Repaired finite-support bounded continuous-kernel normalization for all
  supported Gaussian, Epanechnikov, uniform, and truncated-Gaussian kernels.
  The centered normalization now avoids tail cancellation, uses the platform
  C99 error function where appropriate, and preserves the uniform-density
  limit as the bandwidth tends to infinity.

* Replaced the bounded distribution operator with the mathematically correct
  observation-centered truncated-kernel CDF. Analytic centered interval
  primitives and hoisted invariant work preserve exact support endpoints,
  numerical accuracy, and the finite-support large-bandwidth uniform limit.

* Corrected `transform.bounds = TRUE` initialization for regression initial
  starts and multistarts and for later conditional-distribution multistarts.
  External bandwidth starts are now inverse-mapped before transformed Powell
  search; the public default remains `FALSE`.

* Clarified optimizer summaries so objective-cache lookups, NOMAD point
  lookups, family-native R-`optim` refinement, and total function-evaluation
  accounting are labelled consistently without implying that their
  denominators are interchangeable.

* Clarified bandwidth-object replay documentation across estimator families.
  Stored bandwidth objects retain search and fit metadata but do not silently
  materialize or expand training data; callers must provide data again where
  the documented replay route requires it.

* Standardized the public `nomad.opts` contract across supported bandwidth
  selectors and hardened native NOMAD callback cleanup, nearest-neighbor cache
  release, and native allocation-dimension guards.

* Reworked `npindexbw()` / `npindex()` internals after 0.70-4. Ichimura and
  Klein-Spady single-index objectives now reuse the established `npreg`
  leave-one-out backend where applicable, preserving the public estimator
  contract while materially improving high-dimensional and local-polynomial
  objective evaluation. The single-index formula/materialization routes were
  tightened so formula-selected variables, explicit bandwidth objects, and
  direct formula calls with explicit smoothing parameters reenter consistently.

* Added user documentation for `npindex()` optimizer choice and beta
  interpretation. The documentation now gives practical guidance on when
  derivative-free Nelder-Mead remains a reasonable low-dimensional default and
  when BFGS is useful for higher-dimensional index searches, and it explains
  the relative interpretation of normalized single-index beta coefficients.

* Hardened single-index fit, evaluation, plotting, and summary behavior.
  Bounded continuous-kernel options are now carried consistently through
  objective, fit, evaluation, variance, and bootstrap routes; Klein-Spady
  confusion-matrix output is guarded against out-of-range fitted values; and
  large-bandwidth shortcuts are restricted to kernel/order combinations whose
  constant-weight approximation was validated.

* Repaired single-index plot bootstrap memory usage. Wild-bootstrap
  self-maps now avoid public-facing `O(n^2)` kernel-weight allocations for the
  ordinary large-sample plot route, while preserving the public fitted and
  plotting contracts.

* Repaired generalized-nearest-neighbor local-polynomial derivative ownership
  across `npreghat()` and related public routes. Mixed-degree local-polynomial
  fits now route available derivative components through the correct owner,
  preserve `apply == H %*% y` contracts, and report unavailable derivative
  components consistently rather than silently applying the wrong operator.

* Repaired local-constant derivative ownership in `npreghat()` so degree-zero
  local-polynomial derivative requests use the analytic local-constant
  derivative contract only when that contract is mathematically available.
  Scalar, matrix, and multi-column apply routes were validated separately.

* Improved conditional mixed-degree local-polynomial gradients for
  `npcdens()` / `npcdist()` and public conditional-gradient accessors. The
  partial-availability contract now matches the regression-family policy:
  available components are returned, unavailable components are represented as
  `NA`, and incoherent metadata fails clearly.

* Repaired formula/data reentry contracts across density, distribution,
  conditional density, conditional distribution, single-index, smooth
  coefficient, partially linear, quantile, conditional-mode, copula, and
  significance-test routes. Explicit estimator `data=` now overrides stored
  bandwidth-object data where that public call shape is supported, formula
  `newdata` is validated against fitted RHS variables, and direct formula
  calls with numeric smoothing parameters no longer misroute formula objects
  as native data.

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
  Repairs include `npscoef` iterated backfit behavior, partially-linear fit
  reentry, quantile inversion/clamping contracts, conditional-mode
  probability/tie handling, copula sample reentry, regression tree/large-h
  predicate alignment, and density/distribution formula/native argument
  consistency.

* Plot-bootstrap memory hardening now covers single-index, conditional-mode,
  and partially-linear plot routes. The default plot evaluation grids remain
  linear in the training sample size, and explicit `neval == ntrain` style
  requests remain user-controlled.

* Categorical-gradient and derivative workspace handling were repaired so
  no-gradient and gradient calls can be mixed safely in one session without
  stale native derivative workspace state affecting later calls.

* Native shadow-object and cache lifecycles were hardened. Conditional-density
  shadow pointers, native objective-cache state, and regression large-h /
  large-lambda caches are cleared at the appropriate top-level lifecycle
  boundary so pointer-keyed helper state cannot leak across independent calls
  or datasets.

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
  documentation/demo/benchmark smoke coverage, and installed package proof
  before release-ready claims.

# np 0.70-4

* Hardened proactive C cleanup paths by clearing extended nearest-neighbor
  alias state in the shared estimator cleanup helper, routing selected
  bandwidth-constructor allocation and unsupported-method failures through
  existing cleanup labels, and guarding large `np_kernelsum()` allocation-size
  products before calling the legacy `alloc_vecd(int)` allocator.

* Added `nomad = "auto"` for local-polynomial degree searches. The automatic
  policy uses exhaustive Powell-style degree search for small
  one-dimensional degree lattices where evidence showed it is more reliable
  than heuristic NOMAD restarts, while preserving NOMAD for larger or
  explicitly requested search surfaces.

* Bandwidth and fit summaries now report cumulative search diagnostics more
  clearly. NOMAD cache output distinguishes repeated point lookups avoided
  by NOMAD from true objective computations; Powell summaries expose
  repeated objective lookups avoided by the package-side cache; hybrid
  NOMAD+Powell and exhaustive Powell timing labels are reported consistently.

* Base-graphics plot scaling and legends have been hardened for multi-panel
  displays. Plots now honor active `mfrow`/`cex` behavior more consistently,
  use role-appropriate legend sizes, and draw factor legends with point glyphs
  matching the plotted estimates.

* Fixed-bandwidth bias-corrected bootstrap plot intervals have been reworked
  across the supported plot families. Bias-corrected centers and intervals
  now share a common centering contract for regression, unconditional density
  and distribution, conditional density and distribution, single-index,
  partially linear, and smooth-coefficient plots where supported.

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

* Single-index formula dispatch with explicit bandwidths has been repaired,
  and single-index bias-corrected plot centers, output payloads, and legends
  have been tightened.

* The native `loadNamespace("crs")` call construction used by CRS-backed
  native search routes is now protected, resolving the RCHK protection finding
  without changing the native NOMAD search contract.

# np 0.70-3

* Added `nplsqreg()`/`nplsqregbw()` as a location-scale quantile-regression
  front end with formula/data and bandwidth-object workflows, scalar/vector
  `tau`, prediction, residual extraction, summaries, and plot routes built on
  the shared quantile plotting engine.

* Supported MADS/NOMAD-backed bandwidth-search routes now use the final native
  `crs` NOMAD C API rather than the retired legacy `snomadr()` fallback.
  This covers the promoted regression, density, distribution, conditional
  density, conditional distribution, smooth-coefficient, single-index,
  partially linear, and location-scale quantile search surfaces where those
  routes support native NOMAD/MADS. The runtime dependency on `crs`
  is now declared in `Imports`, while `LinkingTo` remains for the native
  header.

* Native NOMAD routes now preserve progress best-record reporting, expose
  cache/evaluation diagnostics, honor explicit start and option controls, and
  reject unsupported or indeterminate cache-off settings before solver entry.
  Inadmissible GLP degree candidates are guarded before expensive evaluator
  work.

* `npindexbw(..., method = "ichimura", regtype = c("ll", "lp"))` now reuses
  the established local-polynomial regression objective evaluator for
  fixed-degree and NOMAD degree-search routes. Focused sentinel runs preserved
  selected objective payloads while materially reducing runtime for
  local-linear and local-polynomial Ichimura single-index bandwidth searches.

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
  cache under the same switch; NOMAD solver caching and extended-NN distance
  reuse remain separate mechanisms.

* Continuous large-bandwidth shortcut evaluations can now be disabled with
  `options(np.largeh = FALSE)`, and discrete near-upper-bandwidth shortcut
  evaluations can now be disabled with `options(np.largelambda = FALSE)`.
  Both remain enabled by default. These switches are intended for diagnostic
  timing and reproducibility studies that need to separate tree effects from
  large-bandwidth and large-lambda fast paths without changing the canonical
  dense/tree objective machinery.

* Local-polynomial regression cross-validation now uses a leaner hot
  symmetric weighted-sum loop. Fixed-bandwidth `npregbw(..., regtype = "lp",
  bwmethod = "cv.ls")` objective probes show substantially faster
  local-polynomial CV evaluation while preserving objective values to
  numerical precision; adjacent density bandwidth probes preserve their
  objective values as well.

* Shared weighted outer-product accumulation in `npksum()` now uses a guarded
  BLAS `dgemm` route when the operation is dense, non-permuted, and
  memory-bounded. Focused fixed-bandwidth probes preserve objective values to
  numerical precision while substantially accelerating high-basis
  local-polynomial regression and smooth-coefficient objective rows; small and
  scalar routes remain on the established loop path.

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
  density, and regression objective probes preserved objective values exactly
  while reducing repeated scaling work in shared C hot paths.

* Conditional density and conditional distribution least-squares
  cross-validation now use a size-aware row-block policy for local-polynomial
  objective evaluation. The accepted route keeps the bounded-quadrature cap
  unchanged, bounds transient memory by sample size, and preserves objective
  values to numerical precision while materially reducing evaluator overhead
  for fixed-bandwidth CVLS probes.

* Local-polynomial conditional density maximum-likelihood cross-validation now
  uses the same bounded-memory block machinery for fixed and generalized
  nearest-neighbor bandwidths. Focused `npcdensbw(..., bwmethod = "cv.ml",
  regtype = "lp")` probes preserve objective values and selected bandwidths to
  numerical precision while reducing objective and full-search runtime.

* Large-sample categorical-only regression now has a profile-compressed
  execution route controlled by `options(np.categorical.compress = TRUE)`,
  which is enabled by default. This categorical route is independent of
  `options(np.tree)`. For local constant categorical regression, repeated
  predictor profiles are compressed before fitting, prediction/evaluation,
  standard errors, gradients where meaningful, bandwidth search, hat-helper
  use, and plot bootstrap helpers.
  This preserves the established dense-route numerical contract while greatly
  reducing work for large samples with many repeated
  factor/ordered predictor combinations.

* Categorical-only unconditional density routes now use the same
  profile-compression idea when `options(np.categorical.compress = TRUE)` is
  enabled. The fixed-bandwidth fit/evaluation route preserves dense-route
  fitted/evaluation values while avoiding repeated computation over identical
  categorical profiles, and the bandwidth-search route now uses the same
  compressed support representation for all-categorical data. As with other
  flat categorical search surfaces, selected smoothing parameters may drift by
  optimizer-path amounts while preserving the objective scale.

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
  same-process state from leaking across unrelated data sets.

* Formula variables whose names contain dots, such as `y.irr ~ x`, are no
  longer mistaken for the formula wildcard `.` in conditional density and
  conditional distribution bandwidth routes. The conditional-density bandwidth
  formula route also now expands the actual wildcard form `y ~ .` using the
  supplied `data` frame, matching the conditional-distribution route.

# np 0.70-2

* `npqreg()` is now a fully fledged quantile-regression front end. It
  supports the formula/data workflow, internally computes
  `npcdistbw()` bandwidths when a bandwidth object is not supplied,
  accepts scalar or vector `tau`, reuses selected bandwidths for
  additional quantiles in `plot()`, and exposes the usual S3 surface:
  `fitted()`, `predict()`, `predict(..., se.fit=TRUE)`, `se()`,
  `gradients()`, `summary()`, `print()`, `quantile()`, and `plot()`.

* `npqreg()` prediction now honors the standard `newdata` workflow while
  preserving native `exdat` precedence for compatibility with existing
  `np` call surfaces. Formula-based prediction validates that new data
  contain the required right-hand-side variables.

* `npqreg()` plotting has been expanded for vector quantiles,
  level/gradient displays, ordered predictors, user-specified legends,
  and object-fed plotting of additional `tau` values without recomputing
  cross-validation.

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
  optional `rgl` rendering, with asymptotic and bootstrap intervals for
  copula surfaces where defined.

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
  predictors and tensor/additive/Bernstein local-polynomial bases.

* Core and semiparametric S3 prediction paths have been hardened around
  `newdata`, native evaluation-argument precedence, formula RHS
  validation, and `se.fit` handling.

* Front-end/bandwidth argument hygiene has been tightened so
  estimator-only controls such as `proper` are not forwarded into
  bandwidth selectors that do not accept them.

* Documentation has been refreshed for the promoted `npqreg()`,
  `npconmode()`, and `npcopula()` workflows, including the
  local-polynomial NOMAD route, probability/gradient outputs, plot
  controls, and examples that use the streamlined interfaces.

* The pre-release validation suite was expanded with focused hostile
  argument tests, S3 contract tests, installed/tarball proof scripts,
  and cross-package parity checks for the newly promoted estimator
  families.

# np 0.70-1

* The default multistart cap for bandwidth selection now follows
  `min(2, p)` across the core estimator families, replacing the older
  `min(5, p)` cap. This includes automatic LP degree-search calls when
  `search.engine="nomad"` or `"nomad+powell"` and `nmulti` is not
  supplied explicitly.

* The univariate boundary density helper `npuniden.boundary()` now
  defaults to `nmulti=1`.

* The empirical studies supporting this change are documented under
  `benchmarks/validation/`.

* LP-capable front ends now accept `nomad=TRUE` as a documented
  convenience preset for the recommended automatic NOMAD local-polynomial
  route. Missing settings expand to the same long-form LP/NOMAD defaults
  documented in the bandwidth help pages, and regression formula calls
  such as `npreg(y ~ x, nomad = TRUE)` now carry that shortcut through
  the internally computed bandwidth path.
