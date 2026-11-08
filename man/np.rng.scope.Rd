\name{np_rng_scope}
\alias{np_rng_scope}
\title{Temporary Random Number Generator State}
\description{
Some plotting, testing and bandwidth-search operations temporarily seed R's
random number generator and restore the caller's visible state afterward.
An ordinary unlocked \code{.Random.seed} binding with a supported built-in
state is restored exactly. Fresh sessions and workers with a supported built-in
generator also work without an existing seed binding: the operation restores
its original absence and generator kind. Explicit \code{set.seed()} remains
appropriate when reproducible caller draws are wanted.
}
\details{
A seed vector does not contain all Box--Muller or user-supplied generator state.
Operations that reset or replace a temporary state refuse those configurations
rather than claim to restore them. An active, locked or unrecognized
seed binding is also refused. The package never silently chooses a different
generator. To select a different generator, use \code{RNGkind()} explicitly
and then initialize it with \code{set.seed()}.

This restriction concerns temporary restoration, not all direct random draws.
The existing deterministic MPI transport paths that consume only built-in
uniforms can preserve an unused Box--Muller or custom normal generator without
resetting it when the seed binding exists. An absent binding with a hidden-state
generator is refused even for uniform-only transport, because initializing that
binding can itself reset hidden state. Other replay-specific restrictions
continue to apply.

\code{random.seed = NULL}, where passed through to \code{set.seed()}, is not a
request to continue the current stream. A delayed ordinary seed binding may be
resolved once. Arbitrary callbacks that replace or lock the seed binding during
an operation are not a supported restoration contract.
}
\seealso{\code{\link[base]{Random}}, \code{\link[base]{set.seed}}}
\keyword{utilities}
