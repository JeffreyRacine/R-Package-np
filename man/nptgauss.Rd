\name{nptgauss}
\alias{nptgauss}
\title{Truncated Second-order Gaussian Kernels}
\description{
  \code{nptgauss} provides an interface for setting the truncation
  radius of the truncated second-order Gaussian kernel used
  by \pkg{npRmpi}.
}
\usage{
nptgauss(b)
}

\arguments{
  \item{b}{
    Truncation radius of the kernel.
  }
}
\value{
No return value, called for side effects (sets kernel constants in the
\pkg{npRmpi} C backend).
}

\details{
  \code{nptgauss} allows one to set the truncation radius of the truncated Gaussian kernel used by \pkg{npRmpi}, which defaults to 3. It automatically computes the constants describing the truncated gaussian kernel for the user.
  
We define the truncated gaussion kernel on the interval \eqn{[-b,b]} as:
\deqn{K = \frac{\alpha}{\sqrt{2\pi}}\left(e^{-z^2/2} - e^{-b^2/2}\right)}{K = alpha/sqrt(2*pi)*exp(-b^2/2)*(exp(-z^2/2) - exp(-b^2/2))}

The constant \eqn{\alpha} is computed as:
\deqn{\alpha = \left[\int_{-b}^{b} \frac{1}{\sqrt{2\pi}}\left(e^{-z^2/2} - e^{-b^2/2}\right)\right]^{-1}}{1/integrate(1/sqrt(2*pi)*exp(-b^2/2)*(exp(-z^2/2) - exp(-b^2/2)), -b, b)}

Given these definitions, the derivative kernel is simply:
\deqn{K' = (-z)\frac{\alpha}{\sqrt{2\pi}}e^{-z^2/2}}{-z*alpha/(sqrt(2*pi))*exp(-z^2/2)}

The CDF kernel is:
\deqn{G = \frac{\alpha}{2}\mathrm{erf}(z/\sqrt{2}) + \frac{1}{2} - c_0z}{alpha/2*erf(z/sqrt(2))+0.5-c0*z}

The convolution kernel on \eqn{[-2b,0]} has the general form:
\deqn{H_- = a_0\,\mathrm{erf}(z/2 + b) e^{-z^2/4} + a_1z + a_2\,\mathrm{erf}((z+b)/\sqrt{2}) - c_0}{a0*erf(z/2 + b)*exp(-0.25*z^2) + a1*z + a2*erf((z+b)/sqrt(2)) - c0}

and on \eqn{[0,2b]} it is:
\deqn{H_+ = -a_0\,\mathrm{erf}(z/2 - b) e^{-z^2/4} - a_1z - a_2\,\mathrm{erf}((z-b)/\sqrt{2}) - c_0}{-a0*erf(z/2 - b)*exp(-0.25*z^2) - a1*z - a2*erf((z-b)/sqrt(2)) - c0}

where \eqn{a_0}{a0} is determined by the normalisation condition on H,
\eqn{a_2}{a2} is determined by considering the value of the kernel at
\eqn{z = 0} and \eqn{a_1}{a1} is determined by the requirement that \eqn{H = 0} at \eqn{[-2b,2b]}.
}
\author{
    Tristen Hayfield \email{tristen.hayfield@gmail.com}, Jeffrey S. Racine
    \email{racinej@mcmaster.ca}
}

\examples{
\dontrun{
## Not run in checks: excluded to keep MPI examples stable and check times short.
## The default kernel, a gaussian truncated at +- 3
nptgauss(b = 3.0)
} % enddontrun
}
