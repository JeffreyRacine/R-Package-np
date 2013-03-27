\name{nptgauss}
\alias{nptgauss}
\title{Truncated Second-order Gaussian Kernels}
\description{
  \code{nptgauss} provides an interface for setting the parameters
  describing the truncated second-order Gaussian kernel used
  by \pkg{np}.
}
\usage{
nptgauss(b, alpha, c0, a0, a1, a2)
}

\arguments{
  \item{b}{
    Truncation radius of the kernel.
  }

  \item{alpha}{
    The overall normalization of the kernel multiplied by \eqn{\sqrt{2\pi}}{sqrt(2*pi)}.
  }

  \item{c0}{
    A constant containing the value of
    \deqn{\frac{\alpha}{\sqrt{2\pi}}e^{-b^2/2}}{alpha/sqrt(2*pi)*exp(-b^2/2)}.
  }

  \item{a0}{
    A coefficient of the function describing the convolution kernel, see details.
  }

  \item{a1}{
    A coefficient of the function describing the convolution kernel, see details.
  }
  
  \item{a2}{
    A coefficient of the function describing the convolution kernel, see details.
  }

}

\details{
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

}
\author{
    Tristen Hayfield \email{hayfield@mpia.de}, Jeffrey S. Racine
    \email{racinej@mcmaster.ca}
}

\section{Usage Issues}{
  All of the constants for a given truncation radius must be precomputed
  by the user, using a program such as maxima and passed to \code{nptgauss}.
}

\examples{
## The default kernel, a gaussian truncated at +- 3
nptgauss(b = 3.0,
         alpha = 1.030174731161562,
         c0 = 0.004565578246317041,
         a0 = 0.2993759121518507,
         a1 = 2.0844504723243343E-5,
         a2 = -0.004703343342496734)

}