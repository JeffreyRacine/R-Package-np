\name{np_plot_data}
\alias{np_plot_data}
\alias{print.np_plot_data}
\alias{as.data.frame.np_plot_data}
\alias{plot.np_plot_data}
\title{Common Data and Absolute Intervals Returned by Plot Methods}
\description{
\strong{Bias-corrected bounds are experimental (pre-beta research), incomplete,
and not recommended for general-purpose inference.} See the experimental-status
section below.

Plot methods return an \code{np_plot_data} object for \code{output="data"}
or \code{output="plot-data"}. It contains plotted targets and geometry, not
fitted estimator objects. These methods inspect, extract or replay that data
without fitting, bandwidth selection, bootstrap simulation or RNG use.
}
\usage{
\method{print}{np_plot_data}(x, ...)
\method{as.data.frame}{np_plot_data}(x, row.names = NULL, optional = FALSE,
  panel = NULL, target = NULL, band = NULL, ...)
\method{plot}{np_plot_data}(x, y = NULL, ..., panel = NULL, target = NULL,
  band = NULL, renderer = c("base", "rgl"))
}
\arguments{
\item{x}{A plot-data object.}
\item{y}{Unused; must be \code{NULL}.}
\item{panel}{Panel identifier or one-based index. Extraction requires an explicit
selection when there is more than one panel. Replay defaults to all panels.}
\item{target}{Scalar target identifier or one-based index within the selected
panel. Extraction requires a selection when multiple targets exist, for example
several quantiles or derivative components. Replay defaults to all targets.}
\item{band}{An available band name. Extraction requires a selection if several
bands were computed. Replay defaults to all computed bands. This selects saved
endpoints; it never computes a missing band.}
\item{renderer}{Base graphics or the existing interactive \pkg{rgl} surface
renderer. One-dimensional, coefficient, contour and image replay use base graphics.}
\item{row.names, optional}{Standard data-frame method arguments. \code{optional}
is accepted for generic compatibility; column names remain the schema names.}
\item{...}{Graphical overrides for replay, unused by \code{print}. Additional
extraction arguments are rejected.}
}
\details{
All families use the same endpoint convention: \code{lower} and \code{upper}
are absolute values on the target scale. Do not add the estimate or corrected
estimate to them. The corrected estimate need not be the interval midpoint,
particularly for percentile constructions. A standard error is a genuine
available SE, never an interval width divided by a critical value; its reference
is identified in the corresponding inference record. Unavailable values remain
\code{NA} with explicit status and their evaluation rows are retained.

The object has these components:
\describe{
\item{schema_version}{Integer schema identifier, currently \code{1L}.}
\item{family}{Estimator family name.}
\item{request}{Resolved plot controls and, where applicable, explicitly supplied
argument names. These describe the computation; changing them does not recompute it.}
\item{panels}{Named panels. Each contains \code{grid}, \code{variables},
\code{values}, \code{intervals} and \code{display}. Grid columns have stable
identifiers \code{eval_1}, \code{eval_2}, etc., avoiding collisions with target
field names. \code{variables} maps each coordinate to its original name,
role and label, preserving factor/ordered grid columns. Values are keyed by
\code{row_id} and \code{target_id}, with \code{estimate},
\code{bias_corrected}, \code{std_error} and \code{status}.
Intervals add \code{band}, \code{alpha}, \code{lower}, \code{upper},
\code{family_id} and \code{status}. Display holds only retained plotting geometry
and styles, including requested observations, rugs or bootstrap boxplot summaries.}
\item{targets}{Scalar target identities, quantities, labels, derivative variables
and orders, multi-indices, coefficient terms, response levels and quantile
\code{tau} values where relevant. \code{inference_id} links to method metadata.}
\item{inference}{Compact records of the actual construction, law, calibration,
SE reference, projection and available resolution diagnostics. Categorical
contrast metadata records the lower and upper comparison levels for every grid
row, with orientation \code{upper minus lower}. Full bootstrap
draws, covariance workspaces, fitted models and RNG states are not plot results.}
\item{coverage_families}{The finite-grid event and its explicit panel/row/target
members. Pointwise intervals do not imply joint coverage. Simultaneous and
Bonferroni scopes remain those of the selected method; panels, derivatives and
quantiles are not combined unless that method actually calibrated their union.
No coverage between grid points is implied.}
\item{diagnostics}{Additional plot-result diagnostics when available.}
}
Coefficient displays of plus/minus one standard error use band \code{"std_error"},
\code{alpha=NA}, status \code{"display_only"} and no confidence-family claim.
An empirical copula probability-pair display likewise has no estimated target
or interval. Asking for a missing band is an error.

\code{data} and \code{plot-data} return the same numerical payload under matched
controls; only the latter draws the original estimator plot. \code{plot(x)} on
saved plot data replays its geometry. It does not recompute inferential quantities
or acquire missing fitted state. The original estimator methods retain their
existing annotations and graphical controls.

The schema deliberately replaces plot-only \code{d1}/\code{r1}/\code{cd1}
objects, signed error slots, width arrays and inherited estimator methods.
Fitted-model APIs and control-constructor arguments are unchanged. For plot
arguments, use snake_case spellings such as \code{boot_control},
\code{bias_control}, \code{gaussian_draws} and \code{gradient_order}.
}
\value{
\code{as.data.frame} returns the selected grid and values joined to the selected
absolute interval columns. Interval status is named \code{interval_status} in
this flattened table. Without intervals it returns only coordinates and values.
\code{print} invisibly returns \code{x}. Base replay invisibly returns \code{NULL};
interactive surface replay returns the renderer's widget.
}

\section{Experimental status of bias-corrected inference}{
\strong{All bias-corrected confidence intervals and bands are experimental
(pre-beta research), incomplete, and not recommended for general-purpose
inference.} This status remains in effect unless the documentation explicitly
states otherwise for a particular method and scope. It applies wherever
bias-corrected bounds are available: asymptotic and bootstrap constructions,
including legacy methods, pointwise intervals, Bonferroni intervals and
simultaneous bands, for levels and supported derivatives.

These methods form an ongoing research project. Implementation support,
numerical checks and favorable simulation results in selected settings do not
establish generally valid coverage. Coverage can fail outside the studied
settings, and available options do not constitute a general solution to
bias-corrected inference. Method-specific assumptions, restrictions and known
failures must be considered when experimenting with these bounds.
}

\seealso{\code{\link{plot.np}}, \code{\link{np_boot_control}}}
\examples{
\dontrun{
npRmpi.init(nslaves = 1)
x <- seq(-2, 2, length.out = 30)
bw <- npudensbw(dat = data.frame(x = x), bws = .6,
                bandwidth.compute = FALSE)
p <- plot(bw, errors = "asymptotic", band = "pointwise",
          neval = 7, output = "data")
d <- as.data.frame(p)
d[c("eval_1", "estimate", "lower", "upper")]
p$targets
p$panels[[1]]$variables
plot(p)
npRmpi.quit()
} % enddontrun
}
