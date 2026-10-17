\name{np-data-interfaces}
\alias{np-data-interfaces}
\title{Formula Data and Explicit Native Inputs}
\description{
The role of \code{data} depends on the selected interface. It is never a
request to silently ignore a supplied training sample.
}
\details{
Formula calls use \code{data} for model-frame construction. Supported refits
through a bandwidth object retaining a formula also accept replacement
\code{data}; their established omitted-data and \code{data = NULL} behavior
is unchanged.

Native interfaces use explicitly supplied vectors, matrices or data frames.
Where \code{data} is not a native input, even explicitly supplying
\code{data = NULL} is an error. Omit it, supply inputs such as
\code{xdat = d["x"]} and \code{ydat = d$y}, or evaluate the native call inside
\code{with(d, ...)}. Native inputs use ordinary R evaluation, not a data mask.
Existing argument evaluation and method dispatch are preserved. An error in
a supplied expression may therefore occur before the native-data diagnostic;
MPI methods retain their normal error context.
This includes the native IV interfaces: use \code{with(d, npregiv(y, z, w))}
or the supported formula interface.

Some functions genuinely consume a native argument called \code{data}.
In particular, \code{npcopula(bws, data = observations)}, \code{b.star},
\code{npsdeptest} and \code{npsymtest} retain their documented data inputs.
This policy does not change \code{newdata}, native evaluation arguments, or
which input takes precedence during prediction.

Plotting does not generally provide a replacement-data interface. A supplied
\code{data} argument is rejected when the selected plot route does not consume
it, including explicit \code{NULL}. Refit with the intended data and plot the
result, or use supported explicit plotting inputs. The existing partially
linear formula-bandwidth plot route retains its replacement-data support when
its explicit training inputs are omitted.

The diagnostic for unsupported native data is:
\preformatted{`data=` is supported only with the formula interface. For native calls, supply vectors, matrices or data frames directly, or use `with(data, ...)`.}
}
\seealso{\code{\link{npreg}}, \code{\link{npregiv}}, \code{\link{npcopula}}}
\keyword{models}
