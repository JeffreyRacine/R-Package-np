# Private seeded Gaussian plot draws

The private `.np_udens_d8_pilot_seeded()` entry point owns its seed scope. Call it
directly, never inside `.np_with_seed()` under Box-Muller: R's hidden cached normal
has already been destroyed if that outer wrapper called `set.seed()`.

Inversion and other non-Box-Muller kinds retain the existing in-process seeded
provider. Box-Muller uses the requested built-in uniform generator and seed in one
independent `Rscript --vanilla` process. The child loads neither package nor MPI.
It writes normals in chunks of at most 65536 doubles to an invocation-owned file;
the coordinator reads the exact requested blocks and retains the existing MPI
multiplication owner. No replay by seed or worker-local reseeding is used.
The original unseeded provider still consumes the caller's stream as before.

The child has a 600-second execution bound. Missing Rscript, launch failures,
nonzero exit, incomplete receipts/data, nonfinite values and short reads fail;
there is no replacement generator or serial-computation fallback. Normal/factor
files are closed and removed on exit. The bounded vector workspace is distinct
from the independent R process's base memory and total garbage-collected heap.
The stream file uses eight bytes per requested normal and is not a persistent cache.

Pointwise-only and singleton D8 inference launch no sampler. Public plot dispatch,
core estimators and `se=FALSE` remain unchanged. Public D8 activation, Windows
runtime qualification and the generic seed-wrapper/search/bootstrap caller cone
remain separate: this private entry is not a package-wide Box-Muller repair.
Custom uniform generator replay into an independent process has no supported
contract here; the existing unseeded route is unchanged. No runtime np/npRmpi bridge
or new package dependency is introduced.
